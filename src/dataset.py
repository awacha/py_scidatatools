# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:38:52 2011

@author: andris
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cbook import is_numlike
import scipy.optimize
import fitfunction

class DataSetError(StandardError):
    pass

class DataSet(object):
    """ A general purpose dataset class. It has four special fields: x, y,
        dx (sqrt(variance of x)), dy (sqrt(variance of y)), which can be
        accessed by object.fieldname and object['fieldname'] as well. Basic
        arithmetics is implemented with proper error propagation, i.e. "a+b" or
        "a/b" etc. yield the results expected normally. Other fields can also
        be added (but only via the object['newname']=newfield method), which
        will compare the lengths of the new vector and object['x'] before
        the actual creation of the new field. However, be prepared that custom
        datafields may not be inherited by the results of arithmetic operations.
        Trimming and sanitization works, however.
        
        Plotting is also supported via the plot, semilogx, semilogy, loglog and
        errorbar methods. Optional arguments are forwarded to their matplotlib
        equivalents.
        
        Once plotted, the curves can be trimmed by trimzoomed().
        
        Saving to a text file and sanitization of the dataset are also
        supported.
        
        If you subclass this, you may be interested to redefine self._keytrans.
        This dictionary sets up aliases (even many levels deep) to datafields.
        By default, this is {'x':'_x','y':'_y','dy':'_dy','dx':'_dx'}, i.e. the
        fields with underscores in their name (which are the ones actually
        stored in the object) can also be addressed under the names without
        underscores.
        
        Arithmetic operations are defined for the following pairs:
            class DataSet <-> scalar number
            class DataSet <-> class DataSet *
            class DataSet <-> numpy array **
            class DataSet <-> tuple ***
            
       *   only if (max_i(x1_i-x2_i)<min(xtolerance_1,xtolerance_2) is true
       **  only if the array is of the same length as field _x
       *** the tuple should contain at most 4 elements, in the sequence
              (x,y,dy,dx). Elements not present are assumed to be zero.
              
       Datasets can be converted to a structured numpy array by np.array().
       
       Different transformations can also be applied by the set_transform()
           method. These should be subclasses of Dataset.Transform, and only
           affect plotting and "trimzoomed".
    """
    _xtolerance=1e-9
    _dict=None
    _keytrans={'x':'_x','y':'_y','dx':'_dx','dy':'_dy'}
    _plotaxes=None
    _transform=None
    _MCErrorPropSteps=100
    def __init__(self,x,y,dy=None,dx=None,**kwargs):
#    def __init__(self,*args,**kwargs):
#        print "DataSet.__init__:"
#        print "args:",args
#        print "kwargs:", kwargs
        self._dict={}
#        if 'x' in kwargs.keys():
#            x=kwargs['x']
#            del kwargs['x']
#        else:
#            x=args[0]
#        if 'y' in kwargs.keys():
#            y=kwargs['y']
#            del kwargs['y']
#        else:
#            y=args[1]
#        if 'dx' in kwargs.keys():
#            dx=kwargs['dx']
#            del kwargs['dx']
#        else:
#            if len(args)>3:
#                dx=args[3]
#            else:
#                dx=None
#        if 'dy' in kwargs.keys():
#            dy=kwargs['dy']
#            del kwargs['dy']
#        else:
#            if len(args)>2:
#                dy=args[2]
#            else:
#                dy=None
        self['_x']=np.array(x)
        self['_y']=np.array(y)
        if dx is not None:
            self['_dx']=np.array(dx)
        if dy is not None:
            self['_dy']=np.array(dy)
    def __getattr__(self,key):
        # answer to "self.key"-type queries
        
        #This hack is needed to extend the original attributes by the members of
        # _dict, with a twist of self._keytrans.
        
        #Resolution order is as follows.
        # 0) __getattribute__ should have been called before __getattr__, so it
        #    tries to resolve "key" the original way.
        # 1) check self._dict if it contains the key.
        # 2) if not found, try the translation (self._keytrans[key]).
        # 3) if the above two failed, raise AttributeError
        
        #we need to get these two differently, to avoid infinite loops.
#        selfdict=object.__getattribute__(self,'_dict')
#        selfkeytrans=object.__getattribute__(self,'_keytrans')
        # test if "key" is in self._dict...
        try:
            x=self.__getitem__(key)
        except KeyError,key:
            raise AttributeError(key)
        else:
            return x
    def __delattr__(self,key):
        key=self._unalias_keys(key)
        if key in self._dict.keys():
            self.__delitem__(key)
        object.__delattr__(self,key)
    def __delitem__(self,key):
        key=self._unalias_keys(key)
        del self._dict[key]
    def __setattr__(self,key,value):
        #Resolution order: reverse of that in __getattr__
        # 1) check if self._keytrans[key] is in self._dict
        # 2) if that fails, check if "key" is in self._dict
        # 3) if this fails also, try the original way (object.__setattr__)
        if isinstance(self._dict,dict):
            if isinstance(self._keytrans,dict):
                if (key in self._keytrans.keys()) and ((self._keytrans[key]) in (self._dict.keys())):
                    key=self._keytrans[key]
            if key in self._dict.keys():
                #handle these differently! This is the only purpose of this setattr/getattr hack.
                self.__setitem__(key,value)
                return
        return object.__setattr__(self,key,value)
    def __setitem__(self,key,value):
        key=self._unalias_keys(key)
        if key=='_x':
            self._dict.clear()
            self._dict['_x']=np.array(value)
            return
        if '_x' not in self._dict.keys():
            raise ValueError('"_x" should be the first special field to be added!')
        value=np.array(value)
        if len(value)==len(self._dict['_x']):
            self._dict[key]=value
        else:
            raise ValueError('lengths of key "%s" and "x" do not match!'%(key))
    def __getitem__(self,key):
        key1=self._unalias_keys(key)
        if key1 in self._dict.keys():
            return self._dict[key1]
        elif key1 in ['_dx','_dy','y']:
            self._dict[key1]=np.zeros_like(self._dict['_x'])
            return self._dict[key1]
        else:
            raise KeyError(key)
    def __len__(self):
        return len(self._x)
    def keys(self):
        return self._alias_keys(self._dict.keys())
#        ret=[]
#        for k in self._dict.keys():
 #           ret.append(k)
 #           ret.extend([x for x in self._keytrans.keys() if k==self._keytrans[x]])
 #       return ret
    def plot(self,*args,**kwargs):
        if '_plotx' not in self._dict.keys():
            self._do_transform()
        plt.plot(self._plotx,self._ploty,*args,**kwargs)
        self._plotaxes=plt.gca()
    def errorbar(self,*args,**kwargs):
        if '_plotx' not in self._dict.keys():
            self._do_transform()
        if '_plotdy' not in self.keys():
            dy=None
        else:
            dy=self._plotdy
        if '_plotdx' not in self.keys():
            dx=None
        else:
            dx=self._plotdx
        plt.errorbar(self._plotx,self._ploty,dy,dx,*args,**kwargs)
        self._plotaxes=plt.gca()
    def loglog(self,*args,**kwargs):
        if '_plotx' not in self._dict.keys():
            self._do_transform()
        plt.loglog(self._plotx,self._ploty,*args,**kwargs)
        self._plotaxes=plt.gca()
    def semilogx(self,*args,**kwargs):
        if '_plotx' not in self._dict.keys():
            self._do_transform()
        plt.semilogx(self._plotx,self._ploty,*args,**kwargs)
        self._plotaxes=plt.gca()
    def semilogy(self,*args,**kwargs):
        if '_plotx' not in self._dict.keys():
            self._do_transform()
        plt.semilogy(self._plotx,self._ploty,*args,**kwargs)
        self._plotaxes=plt.gca()
    def save(self,filename,cols=None):
        if cols is None:
            cols=['_x','_y','_dy','_dx']
            #add others: all items in self._dict qualify into this category if
            # neither they nor an alias to them (via self._keytrans) exist in
            # cols.
            colsother=[x for x in self._dict.keys() if not self._key_in_set(x,cols)]
            cols.extend(colsother)
            cols=self._alias_keys(cols)
            cols=self._remove_nonexistent_keys(cols)
            cols=[c for c in cols if not c.startswith('_')]
        f=open(filename,'wt')
        f.write('#%s\n'%' '.join(cols))
        lineformat='\t'.join(['%.16g']*len(cols))+'\n'
        for i in range(len(self[cols[0]])):
            line=tuple([self[x][i] for x in cols])
            f.write(lineformat%line)
        f.close()
    def _key_in_set(self,key,set_):
        #if key is directly in the set, return True
        if key in set_:
            return True
        # find keys to which "key" is an alias
        if self._keytrans[key] in set_:
            return True
        # find the aliases of key
        aliases=[x for x in self._keytrans.keys() if self._keytrans[x]==key]
        if len(set(aliases).intersection(set_))>0:
            return True
        return False
    def _alias_keys(self,k):
        if isinstance(k,list):
            return [self._alias_keys(x) for x in k]
        else:
            while (k in self._keytrans.values()):
                k=[x for x in self._keytrans.keys() if self._keytrans[x]==k][0]
            return k
    def _unalias_keys(self,k):
        if isinstance(k,list):
            return [self._unalias_keys(x) for x in k]
        else:
            while(k in self._keytrans.keys()):
                k=self._keytrans[k]
            return k
    def _remove_nonexistent_keys(self,keys):
        return [k for k in keys if self._unalias_keys(k) in self._dict.keys()]
    def _convert_numcompatible(self,c):
        comp={'x':np.zeros_like(self._x),
              'y':np.zeros_like(self._x),
              'dy':np.zeros_like(self._x),
              'dx':np.zeros_like(self._x)}
        if isinstance(c,self.__class__):
            if len(self._x)!=len(c._x):
                return None
            xtol=min(self._xtolerance,c._xtolerance)
            if max(np.abs(self._x-c._x))<xtol:
                try:
                    comp['x']=c._x;
                    comp['y']=c._y;
                    comp['dy']=c._dy;
                    comp['dx']=c._dx;
                except AttributeError:
                    pass
                return comp
        elif isinstance(c,tuple):
            try:
                comp['_x']=np.ndarray([c[0]]*len(self._x))
                comp['_y']=np.ndarray([c[1]]*len(self._x))
                comp['_dy']=np.ndarray([c[2]]*len(self._x))
                comp['_dx']=np.ndarray([c[3]]*len(self._x))
            except IndexError:
                pass
            return comp
        else:
            if is_numlike(c):
                if np.isscalar(c):
                    comp['x']=self._x
                    comp['y']=np.array([c]*len(self._x))
                    # comp['_dy'] and comp['_dx'] have already been set to zero arrays
                    return comp
                elif isinstance(c,np.ndarray) and len(c)==len(self._x):
                    comp['x']=self._x
                    comp['y']=c
                    return comp
                else:
                    return None
            else:
                return None
    def __add__(self,value):
        # self + value operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        return self.__class__((self._x+comp['x'])*0.5,
                              self._y+comp['y'],
                              np.sqrt(self._dy**2+comp['dy']**2),
                              0.5*np.sqrt(self._dx**2+comp['dx']**2))
    def __radd__(self,value):
        # value + self == self + value
        return self + value
    def __iadd__(self,value):
        # self +=value operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        self._dict['_dx']=0.5*np.sqrt(self._dx**2+comp['dx']**2)
        self._dict['_dy']=np.sqrt(self._dy**2+comp['dy']**2)
        self._dict['_x']=0.5*(self._x+comp['x'])
        self._dict['_y']=self._y+comp['y']
        return self
    def __neg__(self):
        return self.__class__(self._x,-self._y,self._dy,self._dx)
    def __sub__(self,value):
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        return self.__class__((self._x+comp['x'])*0.5,
                              self._y-comp['y'],
                              np.sqrt(self._dy**2+comp['dy']**2),
                              0.5*np.sqrt(self._dx**2+comp['dx']**2))
    def __isub__(self,value):
        # self -=value operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        self._dict['_dx']=0.5*np.sqrt(self._dx**2+comp['dx']**2)
        self._dict['_dy']=np.sqrt(self._dy**2+comp['dy']**2)
        self._dict['_x']=0.5*(self._x+comp['x'])
        self._dict['_y']=self._y-comp['y']
        return self
    def __rsub__(self,value):
        #reduce to self.__neg__().__add__(value)
        return (-self)+value
    def __mul__(self,value):
        # self * value operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        return self.__class__((self._x+comp['x'])*0.5,
                              self._y*comp['y'],
                              np.sqrt(comp['y']**2*self._dy**2+comp['dy']**2*self._y**2),
                              0.5*np.sqrt(self._dx**2+comp['dx']**2))
    def __rmul__(self,value):
        #value * self == self * value
        return self * value
    def __imul__(self,value):
        # self *=value operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        self._dict['_dx']=0.5*np.sqrt(self._dx**2+comp['dx']**2)
        self._dict['_dy']=np.sqrt(self._dy**2*comp['y']**2+self._y**2*comp['dy']**2)
        self._dict['_x']=0.5*(self._x+comp['x'])
        self._dict['_y']=self._y*comp['y']
        return self
    def __div__(self,value):
        # self / value operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        return self.__class__((self._x+comp['x'])*0.5,
                              self._y/comp['y'],
                              np.sqrt(1/comp['y']**2*self._dy**2+comp['dy']**2/comp['y']**4*self._y**2),
                              0.5*np.sqrt(self._dx**2+comp['dx']**2))
    def __rdiv__(self,value):
        # value / self operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        return self.__class__((self._x+comp['x'])*0.5,
                              comp['y']/self._y,
                              np.sqrt(comp['y']**2/self._y**4*self._dy**2+comp['dy']**2/self._y**2),
                              0.5*np.sqrt(self._dx**2+comp['dx']**2))
    def __idiv__(self,value):
        # self /=value operation
        comp=self._convert_numcompatible(value)
        if comp is None:
            return NotImplemented
        self._dict['_dx']=0.5*np.sqrt(self._dx**2+comp['dx']**2)
        self._dict['_dy']=np.sqrt(self._dy**2/comp['y']**2+self._y**2/comp['y']**4*comp['dy']**2)
        self._dict['_x']=0.5*(self._x+comp['x'])
        self._dict['_y']=self._y/comp['y']
        return self
    def values(self):
        return self._dict.values()
    def copy(self):
        return self.__class__(**self)
    def trim(self,xmin=-np.inf,xmax=np.inf,inplace=False):
        ind=(self._x<=xmax)&(self._x>=xmin)
        x=self._x[ind]
        y=self._y[ind]
        dy=self._dy[ind]
        dx=self._dx[ind]
        if inplace:
            self._dict['_x']=x
            self._dict['_y']=y
            self._dict['_dx']=dx
            self._dict['_dy']=dy
            return self
        else:
            return self.__class__(x,y,dy,dx);
    def __array__(self,keys=None):
        """Make a structured numpy array from the current dataset.
        """
        if keys==None:
            keys=self.keys()
            values=self.values()
        else:
            keys1=self._unalias_keys(keys)
            values=[self._dict[k] for k in keys1]
        a=np.array(zip(*values),dtype=zip(keys,[np.double]*len(keys)))
        return a
    def sort(self,order='_x'):
        """Sort the current dataset according to 'order' (defaults to '_x').
        """
        order=self._unalias_keys(order)
        keys=self._unalias_keys(self.keys())
        a=self.__array__(keys)
        std=np.sort(a,order=order)
        for k in self._dict.keys():
            self._dict[k]=std[k]
        return self
    def sanitize(self,accordingto=None,thresholdmin=0,thresholdmax=np.inf,function=None):
        """Do an in-place sanitization on this DataSet, i.e. remove nonfinite and out-of-bound elements.
        
        Inputs:
            accordingto: the field, which should be inspected, or a list of them.
                If None, defaults to self.keys(), i.e. all fields.
            thresholdmin: if the inspected field is smaller than this one,
                the line is disregarded.
            thresholdmax: if the inspected field is larger than this one,
                the line is disregarded.
            function: if this is not None, the validity of the dataline is
                decided from the boolean return value of function(value).
                Should accept a list and return a list of booleans. In this
                case, finiteness and the thresholds are NOT checked.
        """
        if not isinstance(accordingto,list):
            accordingto=[accordingto]
        indices=np.array([True]*len(self._x))
        for a in accordingto:
            a=self._unalias_keys(a)
            if hasattr(function,'__call__'):
                indices&=function(self._dict[a])
            else:
                indices&=np.isfinite(self._dict[a])
                indices&=(self._dict[a]>thresholdmin) & (self._dict[a]<thresholdmax)
        for k in self._dict.keys():
            self._dict[k]=self._dict[k][indices]
        return self
    def trimzoomed(self,inplace=False):
        """Trim dataset according to the current zoom on the last plot.
        
        Inputs:
            inplace: True if the current dataset is to be trimmed. If False,
                a new instance of the same class is returned.
                
        Notes:
            This method is useful to reduce the dataset to the currently viewed
                range. I.e. if this dataset was plotted using one of its plot
                methods (e.g. plot(), loglog(), errorbar(), ...) and the graph
                was zoomed in, then calling this function will eliminate the
                off-graph points.
            You will get undefined results if the axis has been deleted since
                the last plot of this instance.
        """
        if self._plotaxes is None:
            raise ValueError('No plot axes corresponds to this dataset!')
        limits=self._plotaxes.axis()
        indices=(self._plotx>=limits[0])&(self._plotx<=limits[1])&(self._ploty>=limits[2])&(self._ploty<=limits[3])
        newdict={}
        for k in self._dict.keys():
            if inplace:
                self._dict[k]=self._dict[k][indices]
            else:
                newdict[self._alias_keys(k)]=self._dict[k][indices]
        if inplace:
            return self
        else:
            x=self.__class__(**newdict)
            x.set_transform(self._transform)
            return x
    def set_transform(self,transform=None):
        self._transform=transform
        try:
            self._do_transform()
        except:
            self._invalidate_transform()
            self._transform=None
    def get_transform(self):
        return self._transform
    def _invalidate_transform(self):
        for key in [k for k in self._dict.keys() if k.startswith('_plot')]:
            if key in self._dict.keys():
                del self._dict[key]
    def _do_transform(self):
        self._invalidate_transform()
        if self._transform is not None:
            d=self._transform(**self)
        else:
            d=self
        for k in d.keys():
            self._dict['_plot%s'%k]=d[k]
    def evalfunction(self,function,*args,**kwargs):
        x=self._dict['_x']
        y=function(x,*args,**kwargs)
        if self._MCErrorPropSteps>1 and '_dx' in self._dict.keys():
            dy=np.zeros_like(x)
            for i in xrange(self.MCErrorPropSteps):
                dy+=(y-function(x+self._dict['_dx']*np.random.randn(len(x))))**2
            dy=np.sqrt(dy)/(self.MCErrorPropSteps-1)
            dx=self._dict['_dx']
        else:
            dy=None
            dx=None
        ret=self.__class__(x,y,dy,dx)
        ret.set_transform(self._transform)
        return ret
    def plotfitted(self,function,params,dparams=None,chi2=None,dof=None,funcinfo=None):
        if funcinfo is not None:
            funcinfo=funcinfo.copy() # make a copy so we can update it.
        else:
            funcinfo={}
            funcinfo['funcname']='Model function'
            funcinfo['paramnames']=['parameter #%d'%(i+1) for i in range(len(params))]
        cfitted=self.evalfunction(function,*params)
        if 'plotmethod' not in funcinfo.keys():
            funcinfo['plotmethod']='plot'
        self.__getattribute__(funcinfo['plotmethod']).__call__(linestyle=' ',marker='.',color='b')
        cfitted.__getattribute__(funcinfo['plotmethod']).__call__(linestyle='-',marker='',color='r')
        logtext=u"Function: %s\n"%funcinfo['funcname']
        logtext=u"Formula: %s\n"%funcinfo['formula']
        logtext+=u"Parameters:\n"
        for i in range(len(params)):
            logtext+=u"    %s : %g +/- %g\n" % (funcinfo['paramnames'][i],
                                                params[i],dparams[i])
        if chi2 is not None:
            logtext+=u"Reduced chi^2: %g\n"%chi2
        if dof is not None:
            logtext+=u"Degrees of freedom: %d\n"%dof
        if 'logtext' in funcinfo.keys():
            for i in funcinfo['logtext']:
                logtext+=u"%s" % i.__call__(self._dict,params,dparams,chi2,dof)
        plt.text(0.95,0.95,logtext,bbox={'facecolor':'white','alpha':0.6,
                                         'edgecolor':'black'},
                 ha='right',va='top',multialignment='left',
                 transform=plt.gca().transAxes)
    def fit(self,function,parinit=None,funcinfo=None,doplot=True,ext_output=False):
        if isinstance(function,fitfunction.FitFunction):
            funcinfo=function.funcinfo()
            if parinit is None:
                parinit=function.init_arguments(self._dict['_x'],self._dict['_y'])
        #get the initial parameters
        if hasattr(parinit,'__call__'):
            params_initial=parinit(self._dict['_x'],self._dict['_y'])
        else:
            params_initial=parinit
        x=self._dict['_x']
        y=self._dict['_y']
        if '_dy' not in self._dict.keys():
            w=np.ones_like(self._dict['_x'])
            weighting='uniform'
        elif (self._dict['_dy']==0).sum()>0:
            idx=self._dict['_dy']==0
            x=x[idx]
            y=y[idx]
            w=np.ones_like(x)
            weighting='uniform, only points where sigma=0'
        else:
            w=1/self._dict['_dy']
            weighting='instrumental'
        def func(p,x,y,w,f=function):
            return (f(x,*(p.tolist()))-y)*w
        p,cov_x,infodict,mesg,ier=scipy.optimize.leastsq(func,params_initial,
                                                         args=(x,y,w),full_output=1)
#        if ier>4 or ier<1:
#            raise DataSetError(mesg)
        chisquare=(infodict['fvec']**2).sum()
        degrees_of_freedom=len(self._dict['_x'])-len(p)
        if cov_x is None:
            pstd=[np.nan]*len(p)
        else:
            pstd=[ np.sqrt(cov_x[i,i]*chisquare/degrees_of_freedom) for i in range(len(p))]
        sserr=np.sum(((function(x,*(p.tolist()))-y)*w)**2)
        sstot=np.sum((y-np.mean(y))**2*w**2)
        r2=1-sserr/sstot
        if funcinfo is not None and doplot:
            self.plotfitted(function,p,pstd,chisquare,degrees_of_freedom,funcinfo)
        if ext_output:
            infodict['Chi2']=chisquare;
            infodict['dof']=degrees_of_freedom
            infodict['mesg']=mesg
            infodict['ier']=ier
            infodict['weighting']=weighting
            infodict['cov_x']=cov_x
            infodict['R2']=r2
            return p,pstd,infodict
        else:            
            return p,pstd
    @staticmethod
    def new_from_file(filename):
        f=np.loadtxt(filename); # this raises IOError if file cannot be loaded.
        N=f.shape[0]
        if N>0:
            x=f[:,0]
            dx=None
            dy=None
            if f.shape[1]==1: # only one column, x will hold the row numbers
                y=f[:,0]
                x=np.arange(f.shape[0],dtype=np.double)
            if f.shape[1]>1:
                y=f[:,1]
            if f.shape[1]>2:
                dy=f[:,2]
            if f.shape[1]>3:
                dx=f[:,3]
            return DataSet(x,y,dy,dx)
        else:
            raise ValueError('File %s does not contain any data points.'%filename)
        
