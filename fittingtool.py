# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 13:28:12 2011

@author: andris
"""

import Tix as Tk
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import os
import B1python

class FitFunction(object):
    argument_info=[]
    name=""
    formula=""
    defaultargs=[]
    def __init__(self):
        pass
    def __call__(self,x,*args,**kwargs):
        raise NotImplementedError
    def init_arguments(self,*args):
        if len(self.defaultargs)==len(self.argument_info):
            return self.defaultargs
        else:
            return [0]*len(self.argument_info)
    def savedefaults(self,args):
        self.defaultargs=args
    def numargs(self):
        return len(self.argument_info)
        
class FFLinear(FitFunction):
    name="Linear"
    argument_info=[('a','slope'),('b','offset')]
    formula="y(x) = a * x + b"
    def __init__(self):
        FitFunction.__init__(self)
    def __call__(self,x,a,b):
        return a*x+b
 
class FFSine(FitFunction):
    name="Sine"
    argument_info=[('a','amplitude'),('omega','circular frequency'),('phi0','phase'),('y0','offset')]
    formula="y(x) = a * sin( omega * x + phi0 ) + y0"
    def __init__(self):
        FitFunction.__init__(self)
    def __call__(self,x,a,omega,phi0,y0):
        return a*np.sin(omega*x+phi0)+y0

class FFGuinier(FitFunction):
    name="Guinier"
    argument_info=[('G','factor'),('Rg','Radius of gyration')]
    formula="y(x) = G * exp(-x^2*Rg^2/3)"
    def __init__(self,qpow=0):
        FitFunction.__init__(self)
        self.qpow=qpow
        if qpow==2:
            self.name="Guinier (thickness)"
            self.formula="y(x) = G * x^2 * exp(-x^2*Rg^2)"
        elif qpow==1:
            self.name="Guinier (cross-section)"
            self.formula="y(x) = G * x * exp(-x^2*Rg^2/2)"
    def __call__(self,x,G,Rg):
        return G*np.power(x,self.qpow)*np.exp(-x*x*Rg*Rg/(3-self.qpow))

class FFPolynomial(FitFunction):
    name="Polynomial"
    argument_info=[]
    formula=""
    def __init__(self,deg=1):
        FitFunction.__init__(self)
        self.deg=deg
        self.name="Polynomial of order %d"%deg
        self.formula="y(x) = sum_(i=0)^(%d) Ai*x^i" % deg
        self.argument_info=self.argument_info[:]
        for i in range(deg+1):
            self.argument_info.append(('A%d'%i,'Coeff. of the x^%d term'%i))
    def __call__(self,x,*coeffs):
        return np.polynomial.polynomial.polyval(x,coeffs)
    

class FFPowerlaw(FitFunction):
    name="Power-law"
    formula="y(x) = A* x^alpha"
    argument_info=[('A','factor'),('alpha','Exponent')]
    def __init__(self,bgorder=None):
        FitFunction.__init__(self)
        self.bgorder=bgorder 
        if bgorder is not None:
            self.name+=" with order #%d background"%bgorder
            self.formula+=" + sum_(i=0)^(%d) Ai*x^i" %bgorder
            self.argument_info=self.argument_info[:]
            for i in range(bgorder+1):
                self.argument_info.append(('A%d'%i,'Coeff. of the x^%d term'%i))
    def __call__(self,x,A,alpha,*bgcoeffs):
        y1=A*np.power(x,alpha)
        if self.bgorder is not None:
            return y1+np.polynomial.polynomial.polyval(x,bgcoeffs)
        return y1


class FittingTool(Tk.Toplevel):
    class FuncSelector(Tk.LabelFrame):
        functions=[FFLinear(),FFPolynomial(2),FFPolynomial(3),FFSine(),FFGuinier(),FFGuinier(1),FFGuinier(2),FFPowerlaw,FFPowerlaw(0),FFPowerlaw(1)]
        argframerow=[]
        selectedfunction=None
        def __init__(self,master,Narguments=8):
            Tk.LabelFrame.__init__(self,master,label='Select a function')
            self.frame.rowconfigure(1,weight=0)
            self.frame.rowconfigure(2,weight=1)
            self.frame.columnconfigure(0,weight=1)
            self.selector=Tk.ComboBox(self.frame,editable=False,command=self.selectfunction,label='Function:')
            self.selector.grid(row=0,column=0,sticky='NSEW')
            for f in self.functions:
                self.selector.insert(Tk.END,f.name)

            lf=Tk.LabelFrame(self.frame,label='Description')
            lf.grid(row=1,column=0,sticky="NEWS")
            lf.frame.rowconfigure(0,weight=1)
            lf.frame.columnconfigure(1,weight=1)
            Tk.Label(lf.frame,text='Formula:',justify=Tk.LEFT,anchor=Tk.W).grid(row=0,column=0,sticky='NSEW')
            self.formulalabel=Tk.Label(lf.frame,text='',justify=Tk.LEFT,anchor=Tk.W)
            self.formulalabel.grid(row=0,column=1,sticky='NSEW')
    
            lf=Tk.LabelFrame(self.frame,label='Arguments')
            lf.grid(row=2,column=0,sticky="NEWS")
            lf.frame.rowconfigure(0,weight=1)
            lf.frame.columnconfigure(0,weight=1)
            argframe=Tk.ScrolledWindow(lf.frame)
            argframe.grid(row=0,column=0,sticky='NSEW')
            argframe=argframe.window
            argframe.columnconfigure(0,weight=1)
            for i in range(Narguments):
                #name value vary? dep
                self.argframerow.append({})
                self.argframerow[-1]['name']=Tk.Label(argframe,text='argument #%d'%i,justify=Tk.LEFT,anchor=Tk.W)
                self.argframerow[-1]['name'].grid(row=i,column=0,sticky='NSW')
                self.argframerow[-1]['value']=Tk.Entry(argframe)
                self.argframerow[-1]['value'].insert(0,'0')
                self.argframerow[-1]['value'].grid(row=i,column=1,sticky='NSEW')
                self.argframerow[-1]['value'].bind('<Key>',self.validateentries)
                self.argframerow[-1]['vary']=Tk.Checkbutton(argframe,text='Vary?')
                self.argframerow[-1]['vary'].grid(row=i,column=2,sticky='NSW')
                argframe.rowconfigure(i,weight=1)
            self.selector.pick(0)
        def validateentries(self,event):
            if not hasattr(event,'calledafter'):
                event.calledafter=True
                self.after_idle(self.validateentries,event)
                return
            x=event.widget.get()
            try:
                float(x)
            except:
                event.widget['bg']='red'
            else:
                event.widget['bg']='white'
        def getargs(self):
            if self.selectedfunction is None:
                return
            N=self.selectedfunction.numargs()
            vals=[]
            for i in range(N):
                vals.append(float(self.argframerow[i]['value'].get()))
            return vals
        def backupargs(self):
            if self.selectedfunction is None:
                return
            self.selectedfunction.savedefaults(self.getargs())
            pass
        def selectfunction(self,funcname):
            func=self.getfunction(funcname)
            self.backupargs()
            initargs=func.init_arguments()
            self.formulalabel['text']=func.formula
            for i in range(len(self.argframerow)):
                if i<len(func.argument_info):
                    for x in self.argframerow[i].values():
                        x['state']='normal'
                    self.argframerow[i]['name']['text']='%s (%s)'%(func.argument_info[i][0],func.argument_info[i][1])
                    self.argframerow[i]['vary'].select()
                    self.argframerow[i]['value'].delete(0,Tk.END)
                    self.argframerow[i]['value'].insert(0,initargs[i])
                else:
                    for x in self.argframerow[i].values():
                        x['state']='disabled'
                    self.argframerow[i]['name']['text']='unused'
                    self.argframerow[i]['vary'].deselect()
                    self.argframerow[i]['value'].delete(0,Tk.END)
            self.selectedfunction=func
        def getfunction(self,funcname=None):
            if funcname is None:
                funcname=self.selector['value']
            return [f for f in self.functions if f.name==funcname][0]
    class DatasetSelector(Tk.LabelFrame):
        def __init__(self,master):
            Tk.LabelFrame.__init__(self,master,label='Dataset')
            self.frame.rowconfigure(0,weight=1)
            self.frame.columnconfigure(0,weight=1)
            b=Tk.Button(self.frame,text='Load file...',command=self.loadfile)
            b.grid(row=0,column=0,columnspan=2,sticky='NSEW')
            Tk.Label(self.frame,text='File:',justify=Tk.LEFT,anchor=Tk.W).grid(row=1,column=0,sticky='NSW')
            Tk.Label(self.frame,text='Directory:',justify=Tk.LEFT,anchor=Tk.W).grid(row=2,column=0,sticky='NSW')
            Tk.Label(self.frame,text='No. of points:',justify=Tk.LEFT,anchor=Tk.W).grid(row=3,column=0,sticky='NSW')
            Tk.Label(self.frame,text='Error bars:',justify=Tk.LEFT,anchor=Tk.W).grid(row=4,column=0,sticky='NSW')
            Tk.Label(self.frame,text='Left:',justify=Tk.LEFT,anchor=Tk.W).grid(row=5,column=0,sticky='NSW')
            Tk.Label(self.frame,text='Right:',justify=Tk.LEFT,anchor=Tk.W).grid(row=6,column=0,sticky='NSW')
            self.filelabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.filelabel.grid(row=1,column=1,sticky='NSEW')
            self.dirlabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.dirlabel.grid(row=2,column=1,sticky='NSEW')
            self.npointslabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.npointslabel.grid(row=3,column=1,sticky='NSEW')
            self.errorbarslabel=Tk.Label(self.frame,text='--',justify=Tk.LEFT,anchor=Tk.W)
            self.errorbarslabel.grid(row=4,column=1,sticky='NSEW')
            self.leftentry=Tk.Entry(self.frame)
            self.leftentry.grid(row=5,column=1,sticky='NSEW')
            self.leftentry['state']='disabled'
            self.rightentry=Tk.Entry(self.frame)
            self.rightentry.grid(row=5,column=1,sticky='NSEW')
            self.rightentry['state']='disabled'
            
        def loadfile(self,name):
            dataset=B1python.readintfile(name)
            if len(dataset)==0:
                return
            self.filelabel['text']=os.path.split(name)[1]
            self.dirlabel['text']=os.path.split(name)[0]
            self.npointslabel['text']='%lu'%len(dataset.q)
            self.errorbarslabel=
            self.leftentry['state']='normal'
            self.leftentry.delete(0,Tk.END)
            self.leftentry.insert(0,'%lf'%(dataset.q.min()))
            self.rightentry['state']='normal'
            self.rightentry.delete(0,Tk.END)
            self.rightentry.insert(0,'%lf'%(dataset.q.max()))
    def __init__(self,root,figure=None):
        Tk.Toplevel.__init__(self,root)
        self.winfo_toplevel().wm_protocol('WM_DELETE_WINDOW',self.quitprogram)
        self.winfo_toplevel().wm_title('Fit tool')
        if figure is None:
            self.figure=plt.figure()
            self.figure.show()
            tkw=self.figure.canvas.get_tk_widget()
            tkw.winfo_toplevel().wm_protocol("WM_DELETE_WINDOW",self.clearandiconizefig)
        self.dss=self.DatasetSelector(self)
        self.dss.grid(sticky='NSEW')
        self.fs=self.FuncSelector(self,Narguments=10)
        self.fs.grid(sticky='NSEW')
        self.rowconfigure(0,weight=1)
        self.columnconfigure(0,weight=1)
    def plot(self,*args,**kwargs):
        self.figure.gca().plot(*args,**kwargs)
        plt.draw()
    def clearfigure(self):
        self.figure.clf()
        plt.draw()
    def clearandiconizefig(self):
        self.clearfigure();
        self.figure.canvas.get_tk_widget().winfo_toplevel().wm_iconify()
        
    def quitprogram(self):
        quit()
        
root=Tk.Tk()
root.withdraw()

mw=FittingTool(root)
mw.mainloop()