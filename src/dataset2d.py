# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 14:11:08 2011

@author: -
"""

import paramstructure
import matplotlib.pyplot as plt
import numpy as np
import scipy.io


class DataSet2DError(Exception):
    pass

class DataSet2D(object):
    """Two-dimensional data and error matrices with optional abscissa and mask.
    
    x is the column coordinate, y is the row (imshow will plot this intuitively).
    
    
    """
    def __init__(self,z,dz='sqrt',x=None,y=None,mask=None,params={}):
        self._xvec=None
        self._yvec=None
        self._dz=None
        self._mask=None
        self.z=np.array(z)
        self.x=x # _setx will be called
        self.y=y # _sety will be called
        if dz is None:
            self._dz=None
        elif dz=='sqrt':
            self._dz=np.sqrt(self.z)
        else:
            self.dz=dz
        self.params=paramstructure.ParamStructure(params)
        self.mask=mask
    def shape(self):
        return self.z.shape
    def copy(self):
        return DataSet2D(z=self.z,dz=self.dz,x=self.x,y=self.y,mask=self.mask,params=self.params)
    def save(self,fname):
        mdict={}
        mdict['z']=self.z
        if self.dz is not None:        
            mdict['dz']=self.dz
        if self._xvec is not None:
            mdict['x']=self._xvec
        if self._yvec is not None:
            mdict['y']=self._yvec
        if self._mask is not None:
            mdict['mask']=self._mask
        if fname.lower().endswith('.mat'):
            scipy.io.savemat(fname,mdict,appendmat=False,oned_as='row')
        else:
            np.savez(fname,**mdict)
    def imshow(self,maskalpha=0.7,*args,**kwargs):
        if not 'origin' in kwargs.keys():
            kwargs['origin']='upper'
        if not 'extent' in kwargs.keys():
            kwargs['extent']=(self.x.min(),self.x.max(),self.y.max(),self.y.min())
        if not 'interpolation' in kwargs.keys():
            kwargs['interpolation']='nearest'
        if self.mask is not None:
            m=self.z.copy()
            m[-np.isfinite(self.z)]=np.min(self.z[np.isfinite(self.z)])
        else:
            m=self.z
        plt.imshow(m,*args,**kwargs)
        if self.mask is not None:
            white=np.ones((self.mask.shape[0],self.mask.shape[1],4))
            white[:,:,3]=(-self.mask).astype('float')*maskalpha
            plt.imshow(white,*args,**kwargs)
        
    def _calcabscissa(self):
        self._xmat,self._ymat=np.meshgrid(self._xvec,self._yvec)
    def getx(self):
        if self._xvec is None:
            self._xvec=np.arange(self.z.shape[1]).astype(np.double)
        return self._xvec
    def setx(self,x):
        if x is None:
            del self.x
            self._xvec=None
            return
        x=np.array(x)
        if x.ndim==1 and len(x)==self.z.shape[1]:
            del self.x
            self._xvec=x
        else:
            raise DataSet2DError('x should be a one-dimensional vector of length M (the number of columns in the matrix)')
    def delx(self):
        if hasattr(self,'_xmat'):
            del self._xmat
        if hasattr(self,'_xvec'):
           del self._xvec
        self._xvec=None
    x=property(getx,setx,delx,'Abscissa (column coordinates)')
    def gety(self):
        if self._yvec is None:
            self._yvec=np.arange(self.z.shape[0]).astype(np.double)
        return self._yvec
    def sety(self,y):
        if y is None:
            del self.y
            self._yvec=None
            return
        y=np.array(y)
        if y.ndim==1 and len(y)==self.z.shape[0]:
            del self.y
            self._yvec=y
        else:
            raise DataSet2DError('y should be a one-dimensional vector of length N (the number of rows in the matrix)')
    def dely(self):
        if hasattr(self,'_ymat'):
            del self._ymat
        if hasattr(self,'_yvec'):
            del self._yvec
        self._yvec=None
    y=property(gety,sety,dely,'Abscissa (row coordinates)')
    def getdz(self):
        return self._dz
    def setdz(self,dz):
        if dz is None:
            self._dz=None
            return
        dz=np.array(dz)
        if dz.shape!=self.z.shape:
            raise DataSet2DError('dz and z should have the same size!')
        self._dz=dz
    def deldz(self):
        if hasattr(self,'_dz'):
            del self._dz
        self._dz=None
    dz=property(getdz,setdz,deldz,'Error matrix')
    def getmask(self):
        return self._mask
    def setmask(self,mask):
        del self.mask
        if mask is None:
            self._mask=None
            return
        mask=np.array(mask).astype('bool')
        if mask.shape!=self.z.shape:
            raise DataSet2DError('mask and z should have the same size!')
        self._mask=mask
    def delmask(self):
        if hasattr(self,'_mask'):
            del self._mask
        self._mask=None
    mask=property(getmask,setmask,delmask,"Mask matrix (bool type, True is non-masked, False is masked)")
    def getxmat(self):
        if not hasattr(self,'_xmat'):
            self._calcabscissa()
        return self._xmat
    def getymat(self):
        if not hasattr(self,'_ymat'):
            self._calcabscissa()
        return self._ymat
    def rescale(self,origin,pixelsize,inplace=False):
        if not inplace:
            obj=self.copy()
        else:
            obj=self
        if np.isscalar(origin):
            origin=[origin,0]
        if np.isscalar(pixelsize):
            pixelsize=[pixelsize]*2
        obj._xvec-=origin[0]
        obj._yvec-=origin[1]
        obj._xvec*=pixelsize[0]
        obj._yvec*=pixelsize[1]
        return obj
    def sum(self,*args,**kwargs):
        if self.mask is not None:
            s=(self.z*self.mask.astype(np.uint8)).sum(*args,**kwargs)
        else:
            s=self.z.sum(*args,**kwargs)
        if self.dz is not None:
            if self.mask is None:
                ds=np.sqrt((self.dz**2).sum(*args,**kwargs))
            else:
                ds=np.sqrt(((self.dz*self.mask.astype(np.uint8))**2).sum(*args,**kwargs))
            return s,ds
        return s
    def maskinvalid(self):
        if self._mask is None:
            self.mask=np.ones_like(self.z).astype(np.bool)
        self.mask&=np.isfinite(self.z)
        if self.dz is not None:
            self.mask&=np.isfinite(self.z)
    def log(self):
        obj=self.copy()
        obj.z=np.log(self.z)
        if self.dz is not None:        
            obj.dz=np.absolute(self.dz/self.z)
        return obj
    def log10(self):
        obj=self.copy()
        obj.z=np.log10(self.z)
        if self.dz is not None:        
            obj.dz=np.absolute(self.dz/self.z)/np.log(10)
        return obj
    
        
           
            
        
        