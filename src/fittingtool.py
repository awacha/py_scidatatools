#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 13:28:12 2011

@author: andris
"""

import Tix as Tk
import matplotlib
import importlib
import sys

# IPython has a wrapper for matplotlib.use to inhibit switching backends.
# It accepts only one argument however, not two, so we have to work it around.
try:
    matplotlib.use('TkAgg',warn=False)
except TypeError:
    pass

import matplotlib.backends
if not matplotlib.backends.backend=='TkAgg':
    raise RuntimeError('Cannot work with other matplotlib backend (currently %s) than TkAgg.'%matplotlib.backends.backend)
    
import matplotlib.pyplot as plt
import os
import numpy as np
import datetime
import time

if __name__ == '__main__':
    from py_scidatatools.dataset import DataSet
    import py_scidatatools.fitfunction as fitfunction
else:
    from dataset import DataSet
    import fitfunction

class FittingtoolStopIteration(Exception):
    pass

class ValidatedEntry(Tk.Entry):
    def __init__(self,*args,**kwargs):
        Tk.Entry.__init__(self,*args,**kwargs)
        self.bind('<Key>',self.validateentry)
        self.goodbackground=self['bg']
    def validateentry(self,event):
        if not hasattr(event,'calledafter'):
            event.calledafter=True
            self.after_idle(self.validateentry,event)
            return
        x=event.widget.get()
        try:
            float(x)
        except:
            event.widget['bg']='red'
        else:
            event.widget['bg']=self.goodbackground

class ValidatedLabelEntry(Tk.LabelEntry):
    def __init__(self,*args,**kwargs):
        Tk.LabelEntry.__init__(self,*args,**kwargs)
        self.entry.bind('<Key>',self.validateentry)
        self.goodbackground=self['bg']
    def validateentry(self,event):
        if not hasattr(event,'calledafter'):
            event.calledafter=True
            self.after_idle(self.validateentry,event)
            return
        x=event.widget.get()
        try:
            float(x)
        except:
            event.widget['bg']='red'
        else:
            event.widget['bg']=self.goodbackground


class FileselectorHelper(object):
    lastdir=None
    def __init__(self,master,command,title='Select file...',filetypes=[('*','All files')],filenamehinter=None):
        self.callbackcommand=command
        self.filenamehinter=filenamehinter
        self._dialog=Tk.ExFileSelectDialog(master.winfo_toplevel(),command=self.callback)
        self._dialog['title']=title
        self._dialog.fsbox['filetypes']=['"%s" "%s"'%(f[0],f[1]) for f in filetypes]
        self._dialog.fsbox.types.pick(0)
        self.called=False
    def __call__(self,*args,**kwargs):
        if not self.called:
            self._dialog.fsbox.dirlist['value']=FileselectorHelper.lastdir
            self.called=True
        types=list(self._dialog.fsbox.types.slistbox.listbox.get(0,Tk.END))
        curval=types.index(self._dialog.fsbox.types['value'])
        self._dialog.fsbox.types.pick(curval)    
        self._dialog.popup()
        if self.filenamehinter is not None:
            fn=self.filenamehinter.__call__()
            fn=self._dialog.fsbox['pattern'].split()[0].replace('*',fn)
            self._dialog.fsbox.file.entry.delete(0,Tk.END)
            self._dialog.fsbox.file.entry.insert(0,fn)
    def callback(self,filename):
        FileselectorHelper.lastdir=self._dialog.fsbox['directory']
        self.callbackcommand.__call__(filename)

class FittingTool(Tk.Toplevel):
    class Xscaleselector(Tk.LabelFrame):
        def __init__(self,master):
            Tk.LabelFrame.__init__(self,master,label='X scale')
            self.frame.columnconfigure(0,weight=1)
            self.usedataset=Tk.Checkbutton(self.frame,text='Get from dataset',command=self.togglecb)
            self.usedataset.grid(row=0,column=0,sticky='NSW')
            self.xmin=ValidatedLabelEntry(self.frame,label='X min')
            self.xmin.grid(row=1,sticky='NSEW')
            self.xmin.entry.insert(0,'0')
            self.xmax=ValidatedLabelEntry(self.frame,label='X max')
            self.xmax.grid(row=2,sticky='NSEW')
            self.xmax.entry.insert(0,'1')
            self.Nx=ValidatedLabelEntry(self.frame,label='Number of points')
            self.Nx.grid(row=3,sticky='NSEW')
            self.Nx.entry.insert(0,'100')
            self.togglecb()
        def togglecb(self):
            if self.usedataset.getvar(self.usedataset['variable'])=='1':
                newstate='disabled'
            else:
                newstate='normal'
            for x in [self.xmin, self.xmax, self.Nx]:
                x['state']=newstate
        def getxscale(self):
            if self.usedataset.getvar(self.usedataset['variable'])=='1':
                return self.winfo_toplevel().getdataset().x
            else:
                return np.linspace(float(self.xmin.entry.get()),float(self.xmax.entry.get()),float(self.Nx.entry.get()))
    class TransformSelector(Tk.LabelFrame):
        transforms=fitfunction.Factory(fitfunction.Transform)
        def __init__(self,master):
            Tk.LabelFrame.__init__(self,master,label='Transformations')
            self.frame.rowconfigure(0,weight=1)
            self.frame.columnconfigure(0,weight=1)
            self.selector=Tk.ScrolledListBox(self.frame,command=self.selecttransform,height=1,width=1)
            self.selector.grid(sticky='NSEW')
            self.selector=self.selector.listbox
            self.selector['exportselection']=0
            tnames=sorted([t.name for t in self.transforms])
            for t in tnames:
                self.selector.insert(Tk.END,t)
            self.selector.select_set(0)
        def selecttransform(self):
            pass
        def normalize_selection(self):
            n=self.selector.curselection()
            if len(n)==1:
                return n[0]
            elif len(n)==0:
                self.selector.select_set(0)
                return 0
            else:
                for i in n[1:]:
                    self.selector.select_clear(i)
                return n[0]
        def gettransform(self):
            n=self.normalize_selection()
            name=self.selector.get(n)
            return [t for t in self.transforms if t.name==name][0]
            
    class FuncSelector(Tk.LabelFrame):
        functions=fitfunction.Factory(fitfunction.FitFunction)
        argframerow=[]
        selectedfunction=None
        def __init__(self,master,Narguments=8):
            self.fully_constructed=False
            Tk.LabelFrame.__init__(self,master,label='Select a function')
            self.frame.rowconfigure(1,weight=0)
            self.frame.rowconfigure(2,weight=1)
            self.frame.columnconfigure(0,weight=1)
            self.selector=Tk.ComboBox(self.frame,editable=False,command=self.selectfunction,label='Function:')
            self.selector.grid(row=0,column=0,sticky='NSEW')
            self.updatefunctions()
            f=Tk.Frame(self.frame)
            f.grid(row=1,column=0,sticky='NSEW')
            Tk.Label(f,text='Formula:',justify=Tk.LEFT,anchor=Tk.W).grid(row=0,column=0,sticky='NSEW')
            self.formulalabel=Tk.Label(f,text='',justify=Tk.LEFT,anchor=Tk.W)
            self.formulalabel.grid(row=0,column=1,sticky='NSEW')
    
            lf=Tk.LabelFrame(self.frame,label='Arguments')
            lf.grid(row=2,column=0,sticky="NEWS")
            lf.frame.rowconfigure(0,weight=1)
            lf.frame.columnconfigure(0,weight=1)
            argframe=Tk.ScrolledWindow(lf.frame,height=90,width=1)
            argframe.grid(row=0,column=0,sticky='NSEW')
            argframe=argframe.window
            argframe.columnconfigure(1,weight=1)
            Tk.Label(argframe,text='Name',justify=Tk.LEFT,anchor=Tk.W).grid(row=0,column=0,sticky='NSW')
            Tk.Label(argframe,text='Value',justify=Tk.LEFT,anchor=Tk.W).grid(row=0,column=1,sticky='NSEW')
            Tk.Label(argframe,text='Scaling',justify=Tk.LEFT,anchor=Tk.W).grid(row=0,column=2,sticky='NSEW')
            Tk.Label(argframe,text='Vary?',justify=Tk.LEFT,anchor=Tk.W).grid(row=0,column=3,sticky='NSW')
            Tk.Label(argframe,text='Error',justify=Tk.LEFT,anchor=Tk.W).grid(row=0,column=4,sticky='NSEW')
            #argframe.columnconfigure(0,weight=1)
            for i in range(Narguments):
                #name value vary? dep
                self.argframerow.append({})
                self.argframerow[-1]['name']=Tk.Label(argframe,text='argument #%d'%i,justify=Tk.LEFT,anchor=Tk.W)
                self.argframerow[-1]['name'].grid(row=i+1,column=0,sticky='NSEW')
                self.argframerow[-1]['name'].tooltip=Tk.Balloon(self)
                self.argframerow[-1]['name'].tooltip.bind_widget(self.argframerow[-1]['name'],balloonmsg='')
                self.argframerow[-1]['name'].tooltip.config_all('bg','yellow')
                self.argframerow[-1]['name'].tooltip['initwait']=100
                self.argframerow[-1]['value']=ValidatedEntry(argframe,width=5)
                self.argframerow[-1]['value'].insert(0,'0')
                self.argframerow[-1]['value'].grid(row=i+1,column=1,sticky='NSEW')
                self.argframerow[-1]['scaling']=ValidatedEntry(argframe,width=5)
                self.argframerow[-1]['scaling'].insert(0,'1')
                self.argframerow[-1]['scaling'].grid(row=i+1,column=2,sticky='NSEW')
                self.argframerow[-1]['vary']=Tk.Checkbutton(argframe,text='Vary?')
                self.argframerow[-1]['vary'].grid(row=i+1,column=3,sticky='NSW')
                self.argframerow[-1]['error']=Tk.Label(argframe,text='--',justify=Tk.LEFT,anchor=Tk.W)
                self.argframerow[-1]['error'].grid(row=i+1,column=4,sticky='NSW')
                argframe.rowconfigure(i,weight=1)
            self.fully_constructed=True
            self.selector.pick(0)
        def getargs(self,getall=True):
            if self.selectedfunction is None:
                return
            N=self.selectedfunction.numargs()
            vals=[]
            for i in range(N):
                if getall or self.argframerow[i]['vary'].getvar(self.argframerow[i]['vary']['variable'])=='1':
                    vals.append(float(self.argframerow[i]['value'].get()))
            return vals
        def getscaling(self,getall=True):
            if self.selectedfunction is None:
                return
            N=self.selectedfunction.numargs()
            vals=[]
            for i in range(N):
                if getall or self.argframerow[i]['vary'].getvar(self.argframerow[i]['vary']['variable'])=='1':
                    vals.append(float(self.argframerow[i]['scaling'].get()))
        def setargs(self,args,stdargs=None):
            self.backupargs()
            self.updateargs(args,stdargs)
        def updatefunctions(self,filename=None):
            if filename is not None:
                dirname,modulename=os.path.split(filename)
                try:
                    sys.path.insert(0,dirname)
                    modulename=os.path.splitext(modulename)[0]
                    importlib.import_module(modulename)
                    sys.path.remove(dirname)
                except:
                    raise
            if hasattr(self,'selector'):
                self.selector.destroy()
            self.selector=Tk.ComboBox(self.frame,editable=False,command=self.selectfunction,label='Function:')
            self.selector.grid(row=0,column=0,sticky='NSEW')
            self.functions=fitfunction.Factory(fitfunction.FitFunction)
            for f in self.functions:
                self.selector.insert(Tk.END,f.name)
            if self.fully_constructed:
                self.selector.pick(0)
            
        def updateargs(self,args,stdargs=None):
            N=self.selectedfunction.numargs()
            j=0
            for i in range(N):
                if (self.argframerow[i]['vary'].getvar(self.argframerow[i]['vary']['variable'])=='1') or len(args)==N:
                    self.argframerow[i]['value'].delete(0,Tk.END)
                    self.argframerow[i]['value'].insert(0,'%g'%args[j])
                    if stdargs is not None:
                        self.argframerow[i]['error']['text']='%g'%stdargs[j]
                    j+=1
                else:
                    self.argframerow[i]['error']['text']='--'
        def backupargs(self):
            if self.selectedfunction is None:
                return
            self.selectedfunction.pushargs(self.getargs())
        def prevargs(self):
            args=self.getargs(getall=True)
            p=self.selectedfunction.popargs()
            while(all([x==y for x,y in zip(args,p)])):
                p=self.selectedfunction.popargs()
            self.selectedfunction.pushargs(p)
            self.updateargs(p)
        def selectfunction(self,funcname):
            func=self.getfunction(funcname,specialize=False)
            self.backupargs()
            initargs=func.init_arguments()
            self.formulalabel['text']=func.formula
            for i in range(len(self.argframerow)):
                if i<len(func.argument_info):
                    for x in self.argframerow[i].values():
                        x['state']='normal'
                    self.argframerow[i]['name']['text']='%s'%(func.argument_info[i][0])
                    self.argframerow[i]['name'].tooltip.unbind_widget(self.argframerow[i]['name'])
                    self.argframerow[i]['name'].tooltip.bind_widget(self.argframerow[i]['name'],balloonmsg=func.argument_info[i][1])
                    self.argframerow[i]['vary'].select()
                    self.argframerow[i]['value'].delete(0,Tk.END)
                    self.argframerow[i]['value'].insert(0,initargs[i])
                else:
                    for x in self.argframerow[i].values():
                        x['state']='disabled'
                    self.argframerow[i]['name']['text']='unused'
                    self.argframerow[i]['name'].tooltip.unbind_widget(self.argframerow[i]['name'])
                    self.argframerow[i]['name'].tooltip.bind_widget(self.argframerow[i]['name'],balloonmsg='')
                    self.argframerow[i]['vary'].deselect()
                    self.argframerow[i]['value'].delete(0,Tk.END)
            self.selectedfunction=func
        def getfunction(self,funcname=None,specialize=True):
            if funcname is None:
                funcname=self.selector['value']
            func=[f for f in self.functions if f.name==funcname][0]
            if specialize:
                spec={}
                for i in range(len(func.argument_info)):
                    if self.argframerow[i]['vary'].getvar(self.argframerow[i]['vary']['variable'])=='0':
                        spec[func.argument_info[i][0]]=float(self.argframerow[i]['value'].get())
                if len(spec.keys())>0:
                    func=func.specialize(**spec)
            return func
        def evalfunction(self,x):
            return self.getfunction()(x,*(self.getargs(getall=False)))            
        
    class FittingControl(Tk.Frame):
        class Messagewindow(Tk.Toplevel):
            niter=0
            starttime=0
            dobreak=False
            def __init__(self,master=None,*args,**kwargs):
                Tk.Toplevel.__init__(self,master,*args,**kwargs)
                self.withdraw()
                self.wm_title('Fitting, please wait...')
                self.wm_protocol('WM_DELETE',self.breakiter)
                self.transient(self.master.winfo_toplevel())
                self.rowconfigure(0,weight=1)
                self.rowconfigure(1,weight=1)
                l=Tk.Label(self,text='Iterations done up to now:',justify=Tk.LEFT,anchor=Tk.W)
                l.grid(sticky='NSEW')
                self.niterlabel=Tk.Label(self,text='0',justify=Tk.LEFT,anchor=Tk.W)
                self.niterlabel.grid(row=0,column=1,sticky='NSEW')
                l=Tk.Label(self,text='Elapsed time:',justify=Tk.LEFT,anchor=Tk.W)
                l.grid(sticky='NSEW')
                self.elapsedtimelabel=Tk.Label(self,text='0',justify=Tk.LEFT,anchor=Tk.W)
                self.elapsedtimelabel.grid(row=1,column=1,sticky='NSEW')
                b=Tk.Button(self,text='Break after current iteration',command=self.breakiter)
                b.grid(columnspan=2)

            def increment(self):
                self.niter+=1;
                self.niterlabel['text']='%d'%self.niter
                self.elapsedtimelabel['text']='%.1f seconds'%(time.time()-self.starttime)
                self.update()
                if self.dobreak:
                    self.withdraw()
                    return False
                return True
            def reset(self):
                self.niter=0
                self.starttime=time.time()
                self.dobreak=False
                self.deiconify()
                #self.master.winfo_toplevel()['state']='disabled'
            def breakiter(self):
                self.dobreak=True
            def done(self):
                #self.master.winfo_toplevel()['state']='normal'
                self.withdraw()
                
            
        def __init__(self,master):
            Tk.Frame.__init__(self,master)
            #self.frame.rowconfigure(0,weight=1)
            self.columnconfigure(0,weight=1)
            self.columnconfigure(1,weight=1)
            self.columnconfigure(2,weight=1)
            b=Tk.Button(self,text="Plot model",command=self.plotmodel)
            b.grid(row=0,column=0,columnspan=1,sticky='NSEW')
            b=Tk.Button(self,text="Plot dataset",command=self.plotdataset)
            b.grid(row=0,column=1,columnspan=1,sticky='NSEW')
            b=Tk.Button(self,text="Clear",command=self.winfo_toplevel().clearfigure)
            b.grid(row=1,column=0,columnspan=1,sticky='NSEW')
            logfiletypes=[('*.log','Log files(*.log)'),('*.txt','Text files'),('*','All files (*)')]            
            savelog=FileselectorHelper(self,self.winfo_toplevel().savelog,'Select a file...',logfiletypes)
            b=Tk.Button(self,text="Save log",command=savelog)
            b.grid(row=1,column=1,sticky='NSEW')
            b=Tk.Button(self,text="Fit",command=self.fit)
            b.grid(row=2,column=0,columnspan=1,sticky='NSEW')
            self.prevargbutton=Tk.Button(self,text='Prev. arg.',command=self.prevargs)
            self.prevargbutton.grid(row=2,column=1,sticky='NSEW')
            self.fittype=Tk.ComboBox(self)
            self.fittype.grid(row=3,column=0,columnspan=2,sticky='NSEW')
            self.fittype.insert(Tk.END,'Normal fit')
            self.fittype['value']='Normal'
            #self.fittype.insert(Tk.END,'Scaled params.')
            self.fittype.insert(Tk.END,'Apparent fit')
            #self.fittype.insert(Tk.END,'Transformed dataset')
            self.fitepsilon=ValidatedLabelEntry(self,label='Epsilon:')
            self.fitepsilon.entry.insert(0,'0.0')
            self.fitepsilon.grid(row=4,column=0,columnspan=2,sticky='NSEW')

            funcfiletypes=[('*.py','Function files (*.py)'),('*','All files (*)')]            
            loadfcn=FileselectorHelper(self,self.winfo_toplevel().loadfitfunction,'Select a file...',funcfiletypes)
            b=Tk.Button(self,text="Load functions",command=loadfcn)
            b.grid(row=5,column=0,columnspan=2,sticky='NSEW')
            self.messagewindow=self.Messagewindow(self)
        def prevargs(self):
            try:
                self.winfo_toplevel().fs.prevargs()
            except IndexError:
                self.prevargbutton.oldbg=self.prevargbutton['bg']
                self.prevargbutton['bg']='red'
                def dummyfunc(butt=self.prevargbutton):
                    butt['bg']=butt.oldbg
                self.prevargbutton.after(1000,dummyfunc)
        def getfittype(self):
            return self.fittype['value'].split()[0]
        def plotmodel(self):
            # get the standard x-scale (real, not transformed values)
            x=self.winfo_toplevel().getxscale()
            # check if we are doing an apparent fit
            if self.getfittype()=='Apparent':
                # transform the x-scale
                x=self.winfo_toplevel().gettransform().do_transform(x,x)['x']
            #evaluate the function in x
            y=self.winfo_toplevel().fs.evalfunction(x)
            # make a dataset
            ds=DataSet(x,y)
            # if not doing an apparent fit, transform the dataset
            if not self.getfittype()=='Apparent':
                ds.set_transform(self.winfo_toplevel().gettransform())
            else:
                ds.set_transform(None)
            self.winfo_toplevel().plot(ds,'r-')
            self.winfo_toplevel().fs.backupargs()
        def plotdataset(self):
            dataset=self.winfo_toplevel().getdataset()
            dataset.set_transform(self.winfo_toplevel().gettransform())
            self.winfo_toplevel().plot(dataset,'b.')
        def fit(self):
            dataset=self.winfo_toplevel().getdataset()
            fittype=self.getfittype()
            if fittype=='Apparent':
                dataset=dataset.apparent()
            filename=self.winfo_toplevel().getfilename()
            self.winfo_toplevel().log("========== %s ==========\n"%datetime.datetime.now().isoformat())
            self.winfo_toplevel().log("Fitting %s\n"%os.path.split(filename)[1])
            self.winfo_toplevel().log("Folder: %s\n"%os.path.split(filename)[0])
            self.winfo_toplevel().log("Function: %s\n"%self.winfo_toplevel().fs.getfunction().formula)
            self.winfo_toplevel().log("x_min: %g\n"%min(dataset.x))
            self.winfo_toplevel().log("x_max: %g\n"%max(dataset.x))
            self.winfo_toplevel().log("Fit type: %s\n"%fittype)
            self.winfo_toplevel().log("Transform: %s\n"%self.winfo_toplevel().gettransform())
            t0=time.time()
            self._currentfitfunction=self.winfo_toplevel().fs.getfunction()
            self.messagewindow.reset()
            try:
                args=self.winfo_toplevel().fs.getargs(getall=False)
                scaling=self.winfo_toplevel().fs.getscaling(getall=False)
                epsfcn=float(self.fitepsilon.get())
                p,pstd,infodict=dataset.fit(self.fitfunction_hacked,args,doplot=False,ext_output=True,diag=scaling,epsfcn=epsfcn)
            except:
                if self.messagewindow.dobreak:
                    self.winfo_toplevel().log('User break.\n')
                else:
                    raise
            else:
                t1=time.time()
                self.winfo_toplevel().log("Finished in %g secs\n"%(t1-t0))
                self.winfo_toplevel().log("Number of function evaluations: %u\n"%(infodict['nfev']))
                self.winfo_toplevel().log("Reduced Chi2: %g\n"%(float(infodict['Chi2'])/infodict['dof']))
                self.winfo_toplevel().log("DoF: %g\n"%infodict['dof'])
                self.winfo_toplevel().log("R^2: %g\n"%infodict['R2'])
                self.winfo_toplevel().log("Variable\tValue\tError\n")
                self.winfo_toplevel().log("---------------------------------------\n")
                for n,v,e in zip(self.winfo_toplevel().fs.getfunction().argument_info,p,pstd):
                    self.winfo_toplevel().log("%s (%s)\t%g\t%g\n"%(n[0],n[1],v,e))
                self.winfo_toplevel().log("\n")
                self.winfo_toplevel().log("Message: (%d) %s\n"%(infodict['ier'],infodict['mesg']))
                self.winfo_toplevel().log("\n")
                self.winfo_toplevel().fs.setargs(p,pstd);
                self.winfo_toplevel().clearfigure()
                self.plotdataset()
                self.plotmodel()
            finally:
                self.messagewindow.done()
        def savelog(self):
            pass
        def fitfunction_hacked(self,*args,**kwargs):
            self.winfo_toplevel().update()
            if not self.messagewindow.increment():
                raise FittingtoolStopIteration
            return self._currentfitfunction(*args,**kwargs)
    class Messagewindow(Tk.LabelFrame):
        def __init__(self,master):
            Tk.LabelFrame.__init__(self,master,label='Messages')
            self.frame.columnconfigure(0,weight=1)
            self.frame.rowconfigure(0,weight=1)
            self.text=Tk.ScrolledText(self.frame,height=90,width=1)
            self.text.grid(sticky='NSEW')
            self.text.text['state']='disabled' #prevent user editing
            self.text.text['width']=30
        def write(self,string):
            self.text.text['state']='normal'
            self.text.text.insert(Tk.END,string)
            self.text.text['state']='disabled'
            self.text.text.see(Tk.END)
        def save(self,filename):
            f=open(filename,'wt')
            f.write(self.text.text.get('0.0',Tk.END))
            f.close()
        
    class DatasetSelector(Tk.LabelFrame):
        def __init__(self,*args,**kwargs):
            if 'label' not in kwargs.keys():
                kwargs['label']='Dataset Selector'
            if 'dataset' in kwargs.keys():
                dataset=kwargs['dataset']
                del kwargs['dataset']
            else:
                dataset=None
            Tk.LabelFrame.__init__(self,*args,**kwargs)
            self.frame.columnconfigure(1,weight=1)
            
            datafiletypes=[('*','All files (*)'),('*.dat','Data files'),('*.txt','Text files')]            
            loadfsh=FileselectorHelper(self,self.loadfile,'Select a file...',datafiletypes)
        
            b=Tk.Button(self.frame,text='Load file...',command=loadfsh)
            b.grid(row=0,column=0,columnspan=2,sticky='NSEW')
            Tk.Label(self.frame,text='File:',justify=Tk.LEFT,anchor=Tk.W).grid(row=1,column=0,sticky='NSW')
            Tk.Label(self.frame,text='Directory:',justify=Tk.LEFT,anchor=Tk.W).grid(row=2,column=0,sticky='NSW')
            Tk.Label(self.frame,text='No. of points:',justify=Tk.LEFT,anchor=Tk.W).grid(row=3,column=0,sticky='NSW')
            Tk.Label(self.frame,text='Error bars:',justify=Tk.LEFT,anchor=Tk.W).grid(row=4,column=0,sticky='NSW')
            self.filelabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.filelabel.grid(row=1,column=1,sticky='NSEW')
            self.dirlabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.dirlabel.grid(row=2,column=1,sticky='NSEW')
            self.npointslabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.npointslabel.grid(row=3,column=1,sticky='NSEW')
            self.errorbarslabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.errorbarslabel.grid(row=4,column=1,sticky='NSEW')
            f=Tk.Frame(self.frame)
            f.columnconfigure(0,weight=1)
            f.columnconfigure(1,weight=1)
            f.grid(row=5,column=0,columnspan=2,sticky='NSEW')
            self.leftentry=Tk.LabelEntry(f,label='Left:')
            self.leftentry.grid(row=0,column=0,sticky='NSEW')
            self.leftentry['state']='disabled'
            self.rightentry=Tk.LabelEntry(f,label='Right:')
            self.rightentry.grid(row=0,column=1,sticky='NSEW')
            self.rightentry['state']='disabled'
            b=Tk.Button(f,text='Fetch',command=self.getrange)
            b.grid(row=0,column=2,sticky='NSEW')
            if dataset is not None:
                self.setdataset(dataset,'<internal>','<internal>')
        def loadfile(self,name):
            try:
                name=os.path.realpath(name)
                dataset=DataSet.new_from_file(name)
                self.setdataset(dataset,os.path.split(name)[1],os.path.split(name)[0])
            except IOError:
                raise
            except ValueError:
                raise

        def setdataset(self,dataset,name='N/A',dirname='N/A'):
            self.dataset=dataset
            self.filename=name
            self.filelabel['text']=name
            self.dirlabel['text']=dirname
            self.npointslabel['text']='%lu'%len(self.dataset.x)
            if self.dataset.dy.sum()>0:
                self.errorbarslabel['text']='Present'
            else:
                self.errorbarslabel['text']='Absent'
            self.leftentry['state']='normal'
            self.leftentry.entry.delete(0,Tk.END)
            self.leftentry.entry.insert(0,'%lf'%(self.dataset.x.min()))
            self.rightentry['state']='normal'
            self.rightentry.entry.delete(0,Tk.END)
            self.rightentry.entry.insert(0,'%lf'%(self.dataset.x.max()))
        def getdataset(self,trim=True):
            try:
                ltrim=float(self.leftentry.entry.get())
            except ValueError:
                ltrim=-np.inf
            try:
                rtrim=float(self.rightentry.entry.get())
            except ValueError:
                rtrim=np.inf
            return self.dataset.trim(ltrim,rtrim,inplace=False)
        def getfilename(self):
            return self.filename
        def getrange(self):
            self.dataset.set_transform(self.winfo_toplevel().gettransform())
            d1=self.dataset.trimzoomed(axes=self.winfo_toplevel().figure.gca())
            self.leftentry.entry.delete(0,Tk.END)
            self.leftentry.entry.insert(0,'%g'%(d1.x.min()))
            self.rightentry.entry.delete(0,Tk.END)
            self.rightentry.entry.insert(0,'%g'%(d1.x.max()))
            del d1
    def __init__(self,*args,**kwargs):
        if 'figure' in kwargs.keys():
            figure=kwargs['figure']
            del kwargs['figure']
        else:
            figure=None
        if 'dataset' in kwargs.keys():
            dataset=kwargs['dataset']
            del kwargs['dataset']
        else:
            dataset=None
        if 'exitonclose' in kwargs.keys():
            self.exitonclose=kwargs['exitonclose']
            del kwargs['exitonclose']
        else:
            self.exitonclose=False
            
        Tk.Toplevel.__init__(self,*args,**kwargs)
        self.winfo_toplevel().wm_protocol('WM_DELETE_WINDOW',self.quitprogram)
        self.winfo_toplevel().wm_title('Fit tool')
        if figure is None:
            self.figure=plt.figure()
            self.figure.show()
            tkw=self.figure.canvas.get_tk_widget()
            tkw.winfo_toplevel().wm_protocol("WM_DELETE_WINDOW",self.clearandiconizefig)
        else:
            self.figure=figure
        
        f=Tk.Frame(self)
        f.grid(row=0,sticky='NSEW')
        f.columnconfigure(0,weight=1)
        if dataset is not None:
            self.dss=self.DatasetSelector(f,dataset=dataset)
        else:
            self.dss=self.DatasetSelector(f)
        self.dss.grid(row=0,columnspan=1,sticky='NSEW')
        self.fc=self.FittingControl(f)
        self.fc.grid(row=0,column=1,columnspan=1,sticky='NSEW')
        f=Tk.Frame(self)
        f.grid(row=1,sticky='NSEW')
        f.columnconfigure(1,weight=1)
        self.xss=self.Xscaleselector(f)
        self.xss.grid(row=0,column=0,sticky='NSEW')
        self.ts=self.TransformSelector(f)
        self.ts.grid(row=0,column=1,sticky='NSEW')
        pw=Tk.PanedWindow(self,orientation='vertical')
        pw.grid(row=2,columnspan=1,sticky='NSEW')
        toppane=pw.add('top',size=200)
        self.fs=self.FuncSelector(toppane,Narguments=20)
        self.fs.pack(expand=1,fill='both')
        bottompane=pw.add('bottom')
        self.mw=self.Messagewindow(bottompane)
        self.mw.pack(expand=1,fill='both')
        self.rowconfigure(2,weight=1)
        self.columnconfigure(0,weight=1)
    def savelog(self,filename):
        self.mw.save(filename)
    def log(self,text):
        self.mw.write(text)
    def plot(self,*args,**kwargs):
        gca=plt.gca()
        plt.axes(self.figure.gca())
        if isinstance(args[0],DataSet):
            args[0].plot(*(args[1:]),**kwargs)
        else:
            plt.plot(*args,**kwargs)
        plt.draw()
        plt.axes(gca)
    def clearfigure(self):
        self.figure.clf()
        plt.draw()
    def clearandiconizefig(self):
        self.clearfigure();
        self.figure.canvas.get_tk_widget().winfo_toplevel().wm_iconify()
    def getdataset(self,*args,**kwargs):
        return self.dss.getdataset(*args,**kwargs)
    def getxscale(self,*args,**kwargs):
        return self.xss.getxscale(*args,**kwargs)
    def quitprogram(self):
        if self.exitonclose:
            quit()
        else:
            self.destroy()
            self.quit()
    def getfilename(self,*args,**kwargs):
        return self.dss.getfilename(*args,**kwargs)
    def gettransform(self,*args,**kwargs):
        return self.ts.gettransform(*args,**kwargs)
    def setdataset(self,dataset,name=''):
        self.dss.setdataset(dataset,name)
    def loadfitfunction(self,*args,**kwargs):
        print "loadfitfunction"
        print "args:",args
        print "kwargs:",kwargs
        self.fs.updatefunctions(*args,**kwargs)
def startfittingtool(dataset=None,exitonclose=False):
    root=Tk.Tk()
    root.withdraw()
    mw=FittingTool(root,dataset=dataset,exitonclose=exitonclose)
    mw.mainloop()

# if called as a script:
if __name__=="__main__":
    startfittingtool(exitonclose=True)
    #cProfile.run('startfittingtool(True)','fittingtool_profile')