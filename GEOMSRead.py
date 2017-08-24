
def loadGEOMSvar(hdfid,varname,varattributes=False,logger=rootlogger):
  """For a filename or a SD opened file in GEOMS format, load the variable varname (and go to SI units)

  
  Input arguments=pyhdf or coda instance,string of variable to load

  Optional key arguments: varattributes=False (returns the variable attributes)"""
  logger=getlogger(logger,'loadGEOMSvar')
  out=(1,array([]))
  if varattributes: out+=({},)
  varwarning=0;closefile=False
  if type(hdfid)==str:
    try: hdfid=coda.open(hdfid)
    except coda.codac.CodacError: logger.error("CODA could not load %s. Abort"%hdfid);return out
    else: closefile=True
  elif type(hdfid)==SD: pass
  else:
    try: coda.get_description(hdfid) 
    except TypeError: logger.error("Mandotory argument %s should be a filename, a pyhdf.SD or a coda file identifier. Abort"%hdffile);return out;
    hdfid=hdfid
  if type(hdfid)==SD:
    try: 
      outvarname=hdfid.select(varname).get()
    except:
      logger.error("%s is not available"%varname)
      return out
    #get attributes
    varattr=hdfid.select(varname).attributes(full=1)
    varunit=float(varattr['VAR_SI_CONVERSION'][0].split(';')[1])
    varoffset=float(varattr['VAR_SI_CONVERSION'][0].split(';')[0])
    varfillvalue=float(varattr['VAR_FILL_VALUE'][0])
    #print varoffset,varunit,outvarname,varfillvalue
    #give warning if values of outvarname is fill_value
  else: 
    try: outvarname=coda.fetch(hdfid,varname.replace('.','_'))
    except coda.codac.CodacError as ce:
      logger.error("%s is not available (%s):\n\t%s"%(varname,repr(ce),'\n\t'.join(coda.get_field_names(hdfid))))
      return out
    varattr=coda.get_attributes(hdfid,varname.replace('.','_'))
    logger.debug('variable %s attributes=%s'%(varname,varattr.__dict__.keys()))
    vardtype=''.join(varattr.VAR_DATA_TYPE)
    if not vardtype=='STRING':
      varhdfdtype=eval(outvarname.dtype.name);logger.debug('Variable numpy dtype=%s'%outvarname.dtype.name)
      dumstr=''.join(varattr.VAR_SI_CONVERSION).split(';')
      varunit=varhdfdtype(dumstr[1]);varoffset=varhdfdtype(dumstr[0])
      try: varfillvalue=array(varattr.VAR_FILL_VALUE).astype(outvarname.dtype)
      except AttributeError: varfillvalue=None;logger.debug('No fill value attr for %s'%varname)
      except Exception as e: logger.error('Something went wrong with the interpretation of the fill value for %s: %s'%(varname,type(varattr.VAR_FILL_VALUE)));raise e
      logger.debug('Variable %s: shape=%s,varunit=%s,offset=%s,fillvalue=%s'%(varname,outvarname.shape,varunit,varoffset,varfillvalue))
    varattr=varattr.__dict__
    [varattr.update({k:''.join(v)}) for k,v in varattr.items() if type(v)==ndarray and type(v[0])==str]
    #except: pass
  if vardtype!='STRING':
    mask=equal(outvarname,varfillvalue)
    if any(mask): 
      logger.warning(varname+" contains fill values");varwarning=1
      getlogger(logger,varname).debug('dtype=%s, fillvalue=%s'%( outvarname.dtype,varfillvalue))
      outvarname=ma.masked_array(outvarname,mask=mask)
      try: outvarname=outvarname.filled(nan)
      except ValueError as e: pass#print varname,outvarname.dtype;raise e
    outvarname=varoffset+outvarname*varunit
  out=(varwarning,outvarname)
  if varattributes: out+=(varattr,);
  if closefile: coda.close(hdfid)
  return out