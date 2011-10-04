# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 14:11:08 2011

@author: -
"""


# TODO: ErrorValue arithmetics

import paramstructure
from attributealias import AliasedArrayAttributes
from utils import ArithmeticBase
import matplotlib
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
import numpy as np
import scipy.io


class DataSet2DError(Exception):
    pass

class DataSet2D(AliasedArrayAttributes,ArithmeticBase):
    """Two-dimensional data and error matrices with optional abscissa and mask.
    
    x is the column coordinate, y is the row (imshow will plot this intuitively).
    """
    _getorigx=lambda p:p['BeamPosX']
    _getorigy=lambda p:p['BeamPosY']
    _getpixelsizex=lambda p:p['PixelSize']
    _getpixelsizey=lambda p:p['PixelSize']
    def __init__(self,params={},**kwargs):
        ArithmeticBase.__init__(self)
        kwargs['normalnames'] = ['z', 'dz','mask']
        AliasedArrayAttributes.__init__(self, **kwargs)
        self.params=paramstructure.ParamStructure(kwargs)
    def copy_into(self,into):
        """Helper function for copy(): make a deep copy,"""
        if not isinstance(into, DataSet2D):
            raise TypeError('copy_into() cannot copy into other types than \
AliasedVectorAttributes or its subclasses')
        AliasedArrayAttributes.copy_into(self, into)
        into.params=self.params.copy()
    def _init_attribute(self,name):
        #this is called by AliasedAttributes.__getattr__ if 'name' is not found.
        
        #if '_shape' is None or an empty list, we cannot do anything, raise an
        # exception
        name=object.__getattribute__(self,'unalias_keys').__call__(name)
        if name=='_mask':
            shape=object.__getattribute__(self, '_shape')
            if not shape: 
                raise NotImplementedError
            object.__getattribute__(self,'addfield')(name, np.ones(shape, dtype=np.bool), False)
        else:
            return AliasedArrayAttributes._init_attribute(self,name)
    def save_matandnpz(self,fname):
        mdict={}
        for k in self.fields():
            mdict[k]=self.getfield(k)
        if fname.lower().endswith('.mat'):
            scipy.io.savemat(fname,mdict,appendmat=False,oned_as='row')
        else:
            np.savez(fname,**mdict)
    def imshow(self,*args,**kwargs):
        my_kwargs={'maskalpha':0.7}
        for k in kwargs:
            if k in my_kwargs.keys():
                my_kwargs[k]=kwargs[k]
                del kwargs[k]
        if not 'origin' in kwargs.keys():
            kwargs['origin']='upper'
        if not 'extent' in kwargs.keys():
            kwargs['extent']=(self.x.min(),self.x.max(),self.y.max(),self.y.min())
        if not 'interpolation' in kwargs.keys():
            kwargs['interpolation']='nearest'
        if len(args)>0 and isinstance(args[0],matplotlib.axes.Axes):
            ax=args[0]
            args=args[1:]
        else:
            ax=plt.gca()
        m=self._z.copy()
        if not self.hasfield('_mask'):
            # set invalid pixels to the smallest valid value
            m[-np.isfinite(self.z)]=np.min(self.z[np.isfinite(self.z)])
        ax.imshow(m,*args,**kwargs)
        if self.hasfield('_mask'):
            white=np.ones((self._mask.shape[0],self._mask.shape[1],4))
            white[:,:,3]=(-self._mask).astype(np.double)*my_kwargs['maskalpha']
            ax.imshow(white,*args,**kwargs)
#    def rescale(self,origin,pixelsize,inplace=False):
#        if not inplace:
#            obj=self.copy()
#        else:
#            obj=self
#        if np.isscalar(origin):
#            origin=[origin,0]
#        if np.isscalar(pixelsize):
#            pixelsize=[pixelsize]*2
#        obj._xvec-=origin[0]
#        obj._yvec-=origin[1]
#        obj._xvec*=pixelsize[0]
#        obj._yvec*=pixelsize[1]
#        return obj
    def sum(self,*args,**kwargs):
        if self.hasfield('_mask'):
            s=(self._z*self._mask.astype(np.uint8)).sum(*args,**kwargs)
        else:
            s=self._z.sum(*args,**kwargs)
        if self.hasfield('_dz'):
            if self.hasfield('_mask'):
                ds=np.sqrt((self._dz**2).sum(*args,**kwargs))
            else:
                ds=np.sqrt(((self._dz*self._mask.astype(np.uint8))**2).sum(*args,**kwargs))
            return s,ds
        return s
    def maskinvalid(self,func=np.isfinite, checkfields=None):
        """Mask invalid pixels
        
        Inputs:
            func: filtering function. Should accept 2D numpy arrays and return
                numpy bool arrays of the same shape.
            checkfields: list of field names (or a single field name in a
                string). These fields will be checked (if they are present).
                Defaults to all defined fields. The mask field is skipped.
                
        Outputs:
            None, the mask matrix is adjusted to mask out all elements where the
            filtering function is False.
        """
        if checkfields is None:
            checkfields=self.fields()
        if isinstance(checkfields,str):
            checkfields=[checkfields]
        for f in checkfields:
            if self.hasfield(f) and self.unalias_keys(f)!='_mask':
                self._mask&=func(self.getfield(f))
    def log(self):
        obj=self.copy()
        obj._z=np.log(self._z)
        if self.hasfield('_dz'):        
            obj._dz=np.absolute(self._dz/self._z)
        return obj
    def log10(self):
        obj=self.copy()
        obj._z=np.log10(self._z)
        if self.hasfield('_dz'):        
            obj._dz=np.absolute(self._dz/self._z)/np.log(10.)
        return obj
    def _iscompatible(self,obj):
        if not isinstance(obj,DataSet2D):
            return NotImplemented
        else:
            if self.shape()!=obj.shape():
                raise DataSet2DError('Datasets are not compatible in shape')
            if np.std(self.x-obj.x)>self._xprecision:
                raise DataSet2DError('x coordinates of datasets are not compatible')
            if np.std(self.y-obj.y)>self._xprecision:
                raise DataSet2DError('y coordinates of datasets are not compatible')
            return obj.copy()
    def __neg__(self):
        obj=self.copy()
        obj._z=-obj._z
        return obj
    def __iadd__(self,rhs):
        rhs1=self._iscompatible(rhs)
        if rhs1 is not NotImplemented:
            self._z+=rhs1._z
            if self.hasattr('_dz') or rhs1.hasattr('_dz'):
                self._dz=np.sqrt(self._dz**2+rhs1._dz**2)
        else:
            self._z+=rhs
    def __imul__(self,rhs):
        rhs1=self._iscompatible(rhs)
        if rhs1 is not NotImplemented:
            if self.hasattr('_dz') or rhs1.hasattr('_dz'):
                self._dz=np.sqrt(rhs1._z**2*self._dz**2+self._z**2*rhs1._dz**2)
            else:
                pass
            self._z*=rhs1._z
        else:
            self._z*=rhs
    def _recip(self):
        obj=self.copy()
        if obj.hasattr('_dz'):
            obj._dz=np.absolute(obj._dz/obj._z**2)
        obj._z=1./self._z
        return obj
    @classmethod
    def load_mat(cls, filename, errorfilename=None, maskfilename=None, *args, **kwargs):
        """Load a 2D dataset from a Matlab(R) MAT file.
        
        Inputs:
            filename: the name of the file (or an open file-like object, which
                can be fed to scipy.io.loadmat()) See Notes!
            errorfilename: the name of the file containing the error matrix (or
                an open file-like object, which can be fed to scipy.io.loadmat()
                ). See Notes!
            maskfilename:  the name of the file containing the mask matrix (or
                an open file-like object, which can be fed to scipy.io.loadmat()
                ). See Notes!
            <other argumens are forwarded to scipy.io.loadmat>
        Output:
            a new instance of this class.
            
        Notes:
            as MAT files can contain several matrices, it is not enough to give
            the file name. The field name can be given following the syntax
            <filename>[:<fieldname>]. If the field name is not given and the mat
            file contains only one variable, it will be loaded. If the file has
            more fields but no field name is given, an exception is raised.
        """
        FIELDSEPARATOR=':'        
        def findfieldname(filename,sep=FIELDSEPARATOR):
            if filename.rfind(FIELDSEPARATOR)==1:
                # we caught the drive name separator in a M$ filename
                fieldname=None
            elif filename.rfind(FIELDSEPARATOR)<0:
                fieldname=None
            else:
                filename,fieldname=filename.rsplit(FIELDSEPARATOR,1)
            m = scipy.io.loadmat(filename,*args,**kwargs)
            if fieldname is None:
                fieldname=[k for k in m.keys() if not (k.startswith('__') and k.endswith('__'))]
                if len(fieldname)>1:
                    raise ValueError('load_mat: File %s contains more than one fields but no fieldname was defined!'%filename)
                if len(fieldname)==0:
                    raise ValueError('load_mat: No fields in file %s.'%filename)
                fieldname=fieldname[0]
            return m[fieldname]
        my_vars={}
        my_vars['_z']=findfieldname(filename)
        if errorfilename is not None:
            my_vars['_dz']=findfieldname(errorfilename)
        if maskfilename is not None:
            my_vars['_mask']=findfilename(maskfilename)
        return cls(**my_vars)
    def get_origin(self):
        return (self._originx(self.params),self._originy(self.params))
    def get_pixelsize(self):
        return 
        raise NotImplementedError
        