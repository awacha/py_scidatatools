# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:08:35 2011

@author: -
"""

import warnings

def _loadfloat(valuestr,fieldname):
    if not isinstance(fieldname,tuple):
        fieldname=tuple([fieldname])
    valuelist=[float(s.strip()) for s in valuestr.split()]
    if len(fieldname)==1 and len(valuelist)>1:
        valuelist=[valuelist]
    return dict(zip(fieldname,valuelist))

def _loadlong(valuestr,fieldname):
    if not isinstance(fieldname,tuple):
        fieldname=tuple([fieldname])
    valuelist=[long(s.strip()) for s in valuestr.split()]
    if len(fieldname)==1 and len(valuelist)>1:
        valuelist=[valuelist]
    return dict(zip(fieldname,valuelist))

def _loadstring(valuestr,fieldname):
    if not isinstance(fieldname,tuple):
        fieldname=tuple([fieldname])
    valuelist=[s.strip() for s in valuestr.split()]
    if len(fieldname)==1 and len(valuelist)>1:
        valuelist=[valuelist]
    return dict(zip(fieldname,valuelist))

def _loadbool(valuestr,fieldname):
    def _tobool(s):
        s=s.lower()
        try:
            return bool(int(s))
        except ValueError:
            pass
        return (s=='true' or s=='yes' or s=='y' or s=='jawohl')
    if not isinstance(fieldname,tuple):
        fieldname=tuple([fieldname])
    valuelist=[_tobool(s.strip()) for s in valuestr.split()]
    if len(fieldname)==1 and len(valuelist)>1:
        valuelist=[valuelist]
    return dict(zip(fieldname,valuelist))

def _formatbool(val):
    if val:
        return 'y'
    else:
        return 'n'

class ParamStructure(dict):
    _logdict=[('FSN','FSN',long),
             ('Sample-to-detector distance (mm)','Dist',float),
             ('Sample thickness (cm)','Thickness',float),
             ('Sample transmission','Transm',float),
             ('Sample position (mm)','PosSample',float),
             ('Temperature','Temperature',float),
             ('Measurement time (sec)','MeasTime',float),
             ('Scattering on 2D detector (photons/sec)','ScatteringFlux',float),
             ('Dark current subtracted (cps)','dclevel',float),
             ('Dark current FSN','FSNdc',long),
             ('Empty beam FSN','FSNempty',long),
             ('Glassy carbon FSN','FSNref1',long),
             ('Glassy carbon thickness (cm)','Thicknessref1',float),
             ('Energy (eV)','Energy',float),
             ('Calibrated energy (eV)','EnergyCalibrated',float),
             ('Calibrated energy','EnergyCalibrated',float),
             ('Beam x y for integration',('BeamPosX','BeamPosY'),float),
             ('Normalisation factor (to absolute units)','NormFactor',float),
             ('Relative error of normalisation factor (percentage)','NormFactorRelativeError',float),
             ('Beam size X Y (mm)',('BeamsizeX','BeamsizeY'),float),
             ('Pixel size of 2D detector (mm)','PixelSize',float),
             ('Primary intensity at monitor (counts/sec)','Monitor',float),
             ('Primary intensity calculated from GC (photons/sec/mm^2)','PrimaryIntensity',float),
             ('Sample rotation around x axis','RotXsample',float),
             ('Sample rotation around y axis','RotYsample',float),
             ('Sample title','Title',str),
             ('Sample name','Title',str),
             ('Injection between Empty beam and sample measurements?','InjectionEB',bool,_formatbool),
             ('Injection between Glassy carbon and sample measurements?','InjectionGC',bool,_formatbool),
             ('FSNs','FSNs',long),
                ]
    def __init__(self,*args,**kwargs):
        dict.__init__(self,*args,**kwargs)
    def __str__(self,lineprefix='',separator=' : '):
        s=''
        for k in self.keys():
            s+='%s%s%s\n'%(lineprefix,k,self[k])
        return s
    def fromstring(self,string,lineprefix='',separators=':='):
        for s in string.splitlines():
            if s.startswith(lineprefix):
                s=s[len(lineprefix):]
            else:
                continue # skip invalid line
            sep=None
            for sep1 in separators:
                if s.count(sep1)==1:
                    if sep is not None:
                        warnings.warn("Malformed line \"%s\": exactly one separator (: or =) is needed."%s)
                        sep=None
                        break
                    else:
                        sep=sep1
                if s.count(sep1)>1:
                    warnings.warn("Malformed line \"%s\": exactly one separator (: or =) is needed."%s)
                    sep=None
                    break
            if sep is None:
                continue # skip current line
            (left,right)=tuple([a.strip() for a in s.split(sep)])
            try:
                idx=[a[0] for a in self._logdict].index(left)
            except ValueError:
                warnings.warn("Ignoring unknown line: %s"%s)
                continue # skip current line
            if self._logdict[idx][2]==bool:
                d=_loadbool(right,self._logdict[idx][1])
            elif self._logdict[idx][2]==str:
                d=_loadstring(right,self._logdict[idx][1])
            elif self._logdict[idx][2]==long:
                d=_loadlong(right,self._logdict[idx][1])
            elif self._logdict[idx][2]==float:
                d=_loadfloat(right,self._logdict[idx][1])
            else:
                raise NotImplementedError("Unknown type: %s"%str(self._logdict[idx][2]))
            self.update(d)
    def tostring(self,lineprefix='',separator=' : '):
        strout=''
        fieldstowrite1=[a for a in self._logdict if a[1] in self.keys()]
        # remove doubles, i.e. same field name but different logfile key
        fieldstowrite=[]
        for i in range(len(fieldstowrite1)):
            if fieldstowrite1[i][1] not in [a[1] for a in fieldstowrite1[0:i]]:
                fieldstowrite.append(fieldstowrite1[i])
        del fieldstowrite1
        for data in fieldstowrite:
            value=self[data[1]]
            if not (isinstance(value,list) or isinstance(value,tuple)):
                value=[value]
            if len(data)==4:
                if isinstance(data[3],str):
                    thisfieldstr=' '.join([data[3]%v for v in value])
                else:
                    thisfieldstr=' '.join([data[3](v) for v in value])
            else:
                thisfieldstr=' '.join([str(v) for v in value])
            strout+='%s%s%s%s\n'%(lineprefix,data[0],separator,thisfieldstr)
        return strout
    def save(self,fname,*args,**kwargs):
        if hasattr(fname,'write'):
            fname.write(self.tostring(*args,**kwargs))
        else:
            f=open(fname,'wt')
            try:
                f.write(self.tostring(*args,**kwargs))
            finally:
                f.close()
    def load(self,fname,*args,**kwargs):
        if hasattr(fname,'read'):
            self.fromstring(fname.read(),*args,**kwargs)
        else:
            f=open(fname,'rt')
            try:
                self.fromstring(f.read(),*args,**kwargs)
            finally:
                f.close()
    def copy(self):
        obj=self.__class__()
        for k in self.keys():
            if hasattr(self[k],'copy'):
                obj[k]=self[k].copy()
            elif issubclass(self[k],list):
                obj[k]=self[k][:]
            else:
                obj[k]=self[k]
        return obj
