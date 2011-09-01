# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 16:18:46 2011

@author: -
"""

from __future__ import division

import Tix as Tk
from math import *
from functools import partial

class Calculator(Tk.Frame):
    def __init__(self,*args,**kwargs):
        Tk.Frame.__init__(self,*args,**kwargs)
        frm=Tk.Frame(self)
        frm.grid(sticky='NSEW')
        frm.columnconfigure(0,weight=1)
        self.cmdline=Tk.Entry(frm)
        self.cmdline.grid(sticky='NSEW')
        self.cmdline.bind('<Return>',self.evaluate)
        self.cmdline.focus_set()
        self.cmdline['justify']=Tk.RIGHT
        Tk.Button(frm,text='=',command=self.evaluate).grid(row=0,column=1)
        frm=Tk.Frame(self)
        frm.grid(sticky='NSEW')
        nkp=NumericKeypad(frm,self.writecmd)
        nkp.grid(row=1,column=0,sticky='NSEW')
        okp=OperatorKeypad(frm,self.writecmd)
        okp.grid(row=1,column=1,sticky='NSEW')
        fkp=FunctionKeypad(frm,self.writecmd)
        fkp.grid(row=1,column=2,sticky='NSEW')
        self._edited=True
    def evaluate(self,event=None):
        a=eval(self.cmdline.get())
        self.cmdline.delete(0,Tk.END)
        self.cmdline.insert(0,a)
        self._edited=False
    def writecmd(self,cmd):
        if not self._edited and cmd.endswith('()'):
            s=cmd[:-1]+self.cmdline.get()+')'
            self.cmdline.delete(0,Tk.END)
            self.cmdline.insert(0,s)
            self._edited=True
        else:
            self.cmdline.insert(Tk.INSERT,cmd)
            if cmd.endswith('()'):
                self.cmdline.icursor(self.cmdline.index(Tk.INSERT)-1)
        self.cmdline.focus_set()
        
class NumericKeypad(Tk.Frame):
    def __init__(self,master,command,*args,**kwargs):
        Tk.Frame.__init__(self,master,*args,**kwargs)
        self.command=command
        self.rowconfigure(0,weight=1)
        self.rowconfigure(1,weight=1)
        self.rowconfigure(2,weight=1)
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=1)
        self.columnconfigure(2,weight=1)
        Tk.Button(self,text='7',command=partial(self.do,7)).grid(row=0,column=0,sticky='NSEW')
        Tk.Button(self,text='8',command=partial(self.do,8)).grid(row=0,column=1,sticky='NSEW')
        Tk.Button(self,text='9',command=partial(self.do,9)).grid(row=0,column=2,sticky='NSEW')
        Tk.Button(self,text='4',command=partial(self.do,4)).grid(row=1,column=0,sticky='NSEW')
        Tk.Button(self,text='5',command=partial(self.do,5)).grid(row=1,column=1,sticky='NSEW')
        Tk.Button(self,text='6',command=partial(self.do,6)).grid(row=1,column=2,sticky='NSEW')
        Tk.Button(self,text='1',command=partial(self.do,1)).grid(row=2,column=0,sticky='NSEW')
        Tk.Button(self,text='2',command=partial(self.do,2)).grid(row=2,column=1,sticky='NSEW')
        Tk.Button(self,text='3',command=partial(self.do,3)).grid(row=2,column=2,sticky='NSEW')
        Tk.Button(self,text='0',command=partial(self.do,0)).grid(row=3,column=0,columnspan=2,sticky='NSEW')
        Tk.Button(self,text='.',command=partial(self.do,'.')).grid(row=3,column=2,sticky='NSEW')
    def do(self,a):
        self.command(str(a))

class OperatorKeypad(Tk.Frame):
    def __init__(self,master,command,*args,**kwargs):
        Tk.Frame.__init__(self,master,*args,**kwargs)
        self.command=command
        self.rowconfigure(0,weight=1)
        self.rowconfigure(1,weight=1)
        self.rowconfigure(2,weight=1)
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=1)
        self.columnconfigure(2,weight=1)
        Tk.Button(self,text='+',command=partial(self.do,'+')).grid(row=0,column=0,sticky='NSEW')
        Tk.Button(self,text='-',command=partial(self.do,'-')).grid(row=1,column=0,sticky='NSEW')
        Tk.Button(self,text='*',command=partial(self.do,'*')).grid(row=2,column=0,sticky='NSEW')
        Tk.Button(self,text='/',command=partial(self.do,'/')).grid(row=3,column=0,sticky='NSEW')
    def do(self,a):
        self.command(str(a))

class FunctionKeypad(Tk.Frame):
    def __init__(self,master,command,*args,**kwargs):
        Tk.Frame.__init__(self,master,*args,**kwargs)
        self.command=command
        self.rowconfigure(0,weight=1)
        self.rowconfigure(1,weight=1)
        self.rowconfigure(2,weight=1)
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=1)
        self.columnconfigure(2,weight=1)
        Tk.Button(self,text='sin',command=partial(self.do,'sin()')).grid(row=0,column=0,sticky='NSEW')
        Tk.Button(self,text='cos',command=partial(self.do,'cos()')).grid(row=1,column=0,sticky='NSEW')
        Tk.Button(self,text='tan',command=partial(self.do,'tan()')).grid(row=2,column=0,sticky='NSEW')
        Tk.Button(self,text='exp',command=partial(self.do,'exp()')).grid(row=3,column=0,sticky='NSEW')
        Tk.Button(self,text='asin',command=partial(self.do,'asin()')).grid(row=0,column=1,sticky='NSEW')
        Tk.Button(self,text='acos',command=partial(self.do,'acos()')).grid(row=1,column=1,sticky='NSEW')
        Tk.Button(self,text='atan',command=partial(self.do,'atan()')).grid(row=2,column=1,sticky='NSEW')
        Tk.Button(self,text='ln',command=partial(self.do,'log()')).grid(row=3,column=1,sticky='NSEW')
        Tk.Button(self,text='sinh',command=partial(self.do,'sinh()')).grid(row=0,column=2,sticky='NSEW')
        Tk.Button(self,text='cosh',command=partial(self.do,'cosh()')).grid(row=1,column=2,sticky='NSEW')
        Tk.Button(self,text='tanh',command=partial(self.do,'tanh()')).grid(row=2,column=2,sticky='NSEW')
        Tk.Button(self,text='log10',command=partial(self.do,'log10()')).grid(row=3,column=2,sticky='NSEW')
        Tk.Button(self,text='asinh',command=partial(self.do,'asinh()')).grid(row=0,column=3,sticky='NSEW')
        Tk.Button(self,text='acosh',command=partial(self.do,'acosh()')).grid(row=1,column=3,sticky='NSEW')
        Tk.Button(self,text='atanh',command=partial(self.do,'atanh()')).grid(row=2,column=3,sticky='NSEW')
        Tk.Button(self,text='sqrt',command=partial(self.do,'sqrt()')).grid(row=3,column=3,sticky='NSEW')
    def do(self,text):
        self.command(text)
    
    