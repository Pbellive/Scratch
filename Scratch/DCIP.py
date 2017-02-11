import numpy as np

def calcK(Ax, Ay, Bx, By, Mx, My, Nx, Ny):
    dx_am = Ax-Mx
    dx_an = Ax-Nx
    dx_bm = Bx-Mx
    dx_bn = Bx-Nx

    dy_am = Ay-My
    dy_an = Ay-Ny
    dy_bm = By-My
    dy_bn = By-Ny

    ram = np.sqrt(dx_am**2 + dy_am**2)
    ran = np.sqrt(dx_an**2 + dy_an**2)
    rbm = np.sqrt(dx_bm**2 + dy_bm**2)
    rbn = np.sqrt(dx_bn**2 + dy_bn**2)

    K = 1/ram - 1/rbm - 1/ran + 1/rbn
    return K

def appRes(Ax, Ay, Bx, By, Mx, My, Nx, Ny, Vp, In=1.):
    
    K = calcK(Ax, Ay, Bx, By, Mx, My, Nx, Ny)
    rho = (2*np.pi*Vp)/(K*In)
    return rho

def writeDCIP(fname, data):
    formatters = 4*[lambda x: '{:.2f}'.format(x)] + 2*[lambda x: '{:.4e}'.format(x)]

    grp = data.groupby(['Ax', 'Ay', 'Bx', 'By'])
    with open(fname, 'w') as f:
        for k in grp.groups.keys():
            txDat = grp.get_group(k)
            out = [*txDat.iloc[0][['Ax', 'Ay', 'Bx', 'By']].values, txDat.shape[0]]
            txLine = '{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:d}\n'.format(*out)
            f.write(txLine)
            txDat.to_string(f, columns=['Mx', 'My', 'Nx', 'Ny', 'Vp', 'VpUncert'], 
                            formatters=formatters,
                            header=False, index=False)
            f.write('\n')



#------------------------------------2D utilities--------------------------------------

def projecttoline(df,y0):
  # Takes input dataFrame containing 2D pdp data and projects electrode locations onto a line.
  # Line coordinates are relative to y0

  # Fit data to line
  m, b, r, p, err = stats.linregress(df.mx,df.my)

  # Translate data such that line of best fit has intercept zero
  x =  df.loc[:,['mx','nx','ax']].values
  y = (df.loc[:,['my','ny','ay']].values-b)
  y = ((x*x + m*x*y)/(x*x + m*m*x*x))*(m*x)
  y = y + b - y0
  df = df.assign(myLine=y[:,0])
  df = df.assign(nyLine=y[:,1])
  df = df.assign(ayLine=y[:,2])

  return df


# Apparent resistivity functions. First define a switcher function that takes
# dataframe and electrode array type string as input. Switcher function dispatches
# on array type
  
def appRes(df,dtype):
  switcher = {
      'dpdp': appResdpdp,
       'pdp': appRespdp,
  }
  print(switcher.keys())
  func = switcher.get(dtype,lambda x:raise_(Exception('Invalid electrode array type')))
  return func(df)

#Pole dipole
def appRespdp(df):
  ram = np.sqrt((df.ax-df.mx)**2 + (df.ay-df.my)**2)
  ran = np.sqrt((df.ax-df.nx)**2 + (df.ay-df.ny)**2)
  G   = (1/ram - 1/ran)/(2*np.pi)
  return df.vpn/G
      
def appResdpdp(df):
  raise_(Exception('dipole dipole appRes not implemented'))
      
def raise_(ex):
  raise ex


#-------------pseudosection stuff---------------------------------------------------------------

def plotAppResPseudosecDepth(dat,lineN,dtype,lb=-9999,ub=-9999,nLevels=10,nTicks=6):

    expX = 0
    expZ = 0

    # horizontal axis
    minX = 0
    maxX = np.amax(dat[['myLine','nyLine','ayLine']].values)
    #ptsX = np.r_[minX, 0.5*(minX+maxX), maxX, minX]
    
    # Get N-spacing from electrode positions
    #N = (dat.ayLine-dat.nyLine)/(dat.nyLine-dat.myLine)
    midpnt = 0.5*(dat.myLine+dat.nyLine)
    N = abs(dat.ayLine - midpnt)
    minZ = -np.max(N)
    maxZ = -np.min(N)
    #ptsZ = np.r_[maxZ, minZ, maxZ, maxZ]
    
    xv = np.linspace(minX, maxX, 1000)
    zv = np.linspace(minZ, 0, 1000)
    
    X, Z = np.meshgrid(xv,zv)

    points = np.c_[midpnt, -N]
    vals = np.log10(dat['Rho'].values)
    #return interp.griddata(points, vals, (X,Z), method='cubic'),N,midpnt
    #RhoApp = 10**(interp.griddata(points, vals, (X,Z), method='cubic'))
    RhoApp = interp.griddata(points, vals, (X,Z), method='linear')
    
    fig = plt.figure(figsize=(15,6))
    ax = fig.gca()
    ax.get_yaxis().set_tick_params(direction='out')
    ax.get_xaxis().set_tick_params(direction='out')

    ax.set_xlim(minX-120, maxX+120)
    ax.set_ylim(minZ-2, 2)

    # Setup contour levels
    lb = vals.min() if lb == -9999 else np.log10(lb)
    ub = vals.max() if ub == -9999 else np.log10(ub)
    print('bounds are',lb,ub)
    levels  = np.linspace(lb,ub,nLevels)
    cbTicks = np.linspace(lb,ub,nTicks)
    cbTickLabels = ["%.2f" % 10**number for number in cbTicks]
    print('levels:',levels)
    print('ticks:', cbTicks)
    print('tickLabels:', cbTickLabels)

    #cf = plt.contourf(xv,zv,Con,ct,extend='both')
    #cf = plt.contourf(xv,zv,Con,ct,extend='both')
    #plt.contour(xv,zv,Con, [1,2,3], colors='k')
    #cf = plt.contourf(xv,zv,RhoApp,norm=colors.LogNorm(vmin=RhoApp.min(),vmax=RhoApp.max()))
    cf = plt.contourf(xv,zv,RhoApp,levels,extend='both')

    plt.ylabel('N-spacing')
    plt.xlabel('Local Easting (m)')

    cb = fig.colorbar(cf,ticks=cbTicks)
    cb.ax.set_ylabel('Apparent Resistivity ($\Omega$ m)')
    
    cb.ax.set_yticklabels(cbTickLabels)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    if dtype == 'PD':
      ax.set_title('Pole-Dipole Resistivity Section - Line #%i'%lineN)
    elif dtype == 'DP':
      ax.set_title('Dipole-Pole Resistivity Section - Line #%i'%lineN)
    else:
        raise Exception('Invalid data type passed to plotAppResPseudosecDepth')

    plt.tight_layout()

    plt.show()
    fname = 'ConSection_%i_%s.pdf'%(lineN,dtype)
    plt.savefig(fname, format='pdf', dpi=300)
    return fig,ax,cf,cb


#------------------------DCIP2D IO------------------------------------------------------------

# Dict of default option values
optDefaults2d = {'cellSize':20.,
        'DCpercentage':0.10,
        'DCfloor':1e-5,
        'IPpercentage':0.075,
        'IPfloor':0.75,
        'DCbackground':0.5,
        'DCchifact':0.9,
        'DCweights':True,
        'IPbackground':25,
        'IPchifact':0.9,
        'IPweights':True,
        }

def setupDirs(lineN):
    # Make base directory
    basedir = 'L%i'%lineN
    try: shutil.rmtree(basedir)
    except: pass
    os.mkdir(basedir)

    os.mkdir(basedir+'/obs/')
    os.mkdir(basedir+'/topo/')
    os.mkdir(basedir+'/mesh/')
    os.mkdir(basedir+'/weight/')

    Inv = ['All', 'DipolePole', 'PoleDipole']
    # Inv = ['All']

    for inv in Inv:
        os.mkdir(basedir+'/'+inv)        
        os.mkdir(basedir+'/'+inv+'/DC')
        os.mkdir(basedir+'/'+inv+'/IP')

def writeDCObs2d(lineN, dat, percentage=optDefaults2d['DCpercentage'], floor=optDefaults2d['DCfloor'], maxN=20):
    # Make observation files
    basedir = 'L%i'%lineN
    Rx = dat[['Rx1','Rx2']].values
    dat.Rx1 = Rx.min(axis=1)
    dat.Rx2 = Rx.max(axis=1)

    #dat = dat[dat.N<=maxN]

    dat['uncertVP'] = percentage*np.abs(dat.vp) + floor

    #sign = np.ones(dat.shape[0])
    #sign[np.where(dat.DP)] = -1.
    #dat.vp = sign*np.abs(dat.vp)

    #with open(basedir+'/obs/L%i_DC_PoleDipole.obs'%lineN, 'w') as f:
        #f.write('\n')
        #dat.sort_values('Tx2').to_csv(f, 
                     #columns=['Tx2', 'Tx2', 'Rx1', 'Rx2', 'vp', 'uncertVP'],
                     #index=False,
                     #header=False, 
                     #sep=' ')

    with open(basedir+'/obs/L%i_DC_PoleDipole.obs'%lineN, 'w') as f:
        f.write('\n')
        dat[dat.PD].sort('Tx2').to_csv(f, 
                     columns=['Tx2', 'Tx2', 'Rx1', 'Rx2', 'vp', 'uncertVP'],
                     index=False,
                     header=False, 
                     sep=' ')

    with open(basedir+'/obs/L%i_DC_DipolePole.obs'%lineN, 'w') as f:
        f.write('\n')
        dat[dat.DP].to_csv(f, 
                     columns=['Tx2', 'Tx2', 'Rx1', 'Rx2', 'vp', 'uncertVP'],
                     index=False,
                     header=False, 
                     sep=' ')

    with open(basedir+'/obs/L%i_DC_All.obs'%lineN, 'w') as f:
        f.write('\n')
        dat.to_csv(f, 
                     columns=['Tx2', 'Tx2', 'Rx1', 'Rx2', 'vp', 'uncertVP'],
                     index=False,
                     header=False, 
                     sep=' ')

def writeIPObs2d(lineN, dat, percentage=optDefaults2d['IPpercentage'], floor=optDefaults2d['IPfloor'], maxN=20, maxX=np.inf):
    # Make observation files
    basedir = 'L%i'%lineN
    Rx = dat[['Rx1','Rx2']].values
    dat.Rx1 = Rx.min(axis=1)
    dat.Rx2 = Rx.max(axis=1)

    dat = dat[dat.N<=maxN]
    dat = dat[(dat.Rx1<maxX) & (dat.Rx2<maxX) & (dat.Tx2<maxX)]

    dat['uncertChg'] = percentage*np.abs(dat.chg) + floor

    with open(basedir+'/obs/L%i_IP_PoleDipole.obs'%lineN, 'w') as f:
        f.write('\n')
        dat[dat.PD].to_csv(f, 
                     columns=['Tx2', 'Tx2', 'Rx1', 'Rx2', 'chg', 'uncertChg'],
                     index=False,
                     header=False, 
                     sep=' ')

    with open(basedir+'/obs/L%i_IP_DipolePole.obs'%lineN, 'w') as f:
        f.write('\n')
        dat[dat.DP].to_csv(f, 
                     columns=['Tx2', 'Tx2', 'Rx1', 'Rx2', 'chg', 'uncertChg'],
                     index=False,
                     header=False, 
                     sep=' ')

    with open(basedir+'/obs/L%i_IP_All.obs'%lineN, 'w') as f:
        f.write('\n')
        dat.to_csv(f, 
                     columns=['Tx2', 'Tx2', 'Rx1', 'Rx2', 'chg', 'uncertChg'],
                     index=False,
                     header=False, 
                     sep=' ')

def makeTopo2d(lineN):
    basedir = 'L%i'%lineN
    # Make the topo
    topo = pd.read_csv('../topo/lawlessTopo.txt', 
                       skiprows=1, 
                       delim_whitespace=True, 
                       names=['x','y','z'])

    el = pd.read_csv('../raw/june82015/Lawless.csv')
    el['Station'] = el.Label.str.split('_').str[1].apply(float)
    el['Line'] = el.Label.str.split('_').str[0].str[1:-1].apply(float)
    el = el[['Line', 'Station', 'Easting', 'Northing']]
    Fz = LinearNDInterpolator(topo[['x','y']], topo.z)
    line = el[el.Line==lineN]

    St = np.arange(line.Station.min()-500, line.Station.max()+500, 10)

    mEasting = (line.iloc[-1] - line.iloc[0]).Easting/(line.iloc[-1] - line.iloc[0]).Station
    bEasting = line.iloc[0].Easting - mEasting*line.iloc[0].Station

    mNorthing = (line.iloc[-1] - line.iloc[0]).Northing/(line.iloc[-1] - line.iloc[0]).Station
    bNorthing = line.iloc[0].Northing - mNorthing*line.iloc[0].Station


    StEasting = mEasting*St + bEasting
    StNorthing = mNorthing*St + bNorthing


    StZ = Fz(StEasting, StNorthing)

    out = ''
    out += '%i\t%.1f\n' %(StZ.shape[0], StZ.max())
    for i, r in enumerate(np.c_[St, StZ]):
        out += '%.1f\t%.1f\n' %tuple(r)

    with open(basedir+'/topo/topo%i.txt'%lineN, 'w') as f:
        f.write(out)

def makeMesh2d(lineN, dat, cellSize=optDefaults2d['cellSize']):
    basedir = 'L%i'%lineN
    locs = dat[['Tx2', 'Rx1', 'Rx2']].values.flatten()

    nPadd = 15
    expFact = 1.2
    baseCS = cellSize
    corePadd = 200.
    padd = baseCS*1.2**np.r_[0:nPadd]
    coreWidth = locs.max() - locs.min() + 2*corePadd
    nCore = coreWidth/baseCS
    nCoreZ = 1000./baseCS

    offSet = sum(padd) + corePadd - dat.Tx2.min()
    N = nCore + 2*nPadd 
    out = ''
    out += '%i\n' %(nPadd*2 + 1)
    out += '%.2f ' %(-offSet)
    for xi in np.r_[padd[::-1].cumsum()] - offSet:
        out += '%.2f 1\n' %xi
    out += '%.2f %i\n' %(padd.sum() + coreWidth -offSet, nCore)

    for xi in padd.sum() + coreWidth + padd.cumsum() - offSet:
        out += '%.2f 1\n' %xi

    out += '%i\n' %(nPadd + 1)
    out += '0.00 %.2f %i\n' %(nCoreZ*baseCS, nCoreZ) 
    for zi in nCoreZ*baseCS + padd.cumsum():
        out += '%.2f 1\n' %zi

    with open(basedir+'/mesh/mesh%i.txt'%lineN, 'w') as f: f.write(out)

def makeWeights2d(lineN):
    basedir = 'L%i'%lineN
    # Make weights
    wgtStr =  ''
    wgtStr += '../mesh/mesh%i.txt    ! input mesh file\n'%lineN
    wgtStr += '../topo/topo%i.txt    ! topography\n'%lineN
    wgtStr += 'NO_MODEL        ! input model file\n'
    wgtStr += 'USE_LOG\n'
    wgtStr += 'weight.txt      ! output weight file\n'
    wgtStr += '0.01  0.01 ! gradient tolerances in X and Z\n'
    wgtStr += '0.01       ! edge weights\n'
    # wgtStr += '25          ! # of top layers for the X interface weightings\n'
    # wgtStr += '20.0 16.2 13.1 10.7 8.78 7.22 5.98 4.98 4.18 3.55 3.04 2.63 2.30 2.04 1.83 1.66 1.53 1.42 1.34 1.27 1.21 1.17 1.14 1.11 1.08  ! weight values for the top X interfaces\n'
    wgtStr += '8          ! # of top layers for the X interface weightings\n'
    wgtStr += '11 6 3.5 2.25 1.63 1.31 1.16 1.08 1.04  ! weight values for the top X interfaces\n'
    wgtStr += '4          ! # of top layers for the Z interface weightings\n'
    wgtStr += '10 5 2.5 1.25       ! weight values for the top Z interfaces\n'

    with open(basedir+'/weight/find_edges.inp', 'w') as f: f.write(wgtStr)

    os.chdir(basedir+'/weight/')
    call(['/home/patrick/src/dcip2d/find_edges/find_edges','find_edges.inp'])
    os.chdir('../..')

def writeDCInp2d(lineN, background=optDefaults2d['DCbackground'], chifact=optDefaults2d['DCchifact'], useWeights=optDefaults2d['DCweights']):
    basedir = 'L%i'%lineN
    # Make inversion directories
    Inv = ['All', 'DipolePole', 'PoleDipole']
    # Inv = ['All']
    for inv in Inv:
        inpStr = ''
        inpStr += 'MESH FILE ../../mesh/mesh%i.txt\n'%lineN
        inpStr += 'OBS LOC_X ../../obs/L%i_DC_%s.obs\n'%(lineN, inv)
        inpStr += 'CHIFACT %.2f\n'%chifact
        inpStr += 'TOPO FILE ../../topo/topo%i.txt\n'%lineN
        inpStr += 'REF_MOD VALUE %.2e\n'%background
        inpStr += 'WAVE 2.5e-4 1.0 13\n'
        inpStr += 'STORE_ALL_MODELS TRUE\n'
        inpStr += 'USE_MREF TRUE\n'
        inpStr += 'BOUNDS NONE\n'
        inpStr += 'INVMODE CG\n'
        if useWeights:
            inpStr += 'WEIGHT FILES ../../weight/S_weight.txt ../../weight/X_weight.txt ../../weight/Z_weight.txt\n'

        with open(basedir+'/'+inv+'/DC/dcinv2d.inp', 'w') as f: f.write(inpStr)

def writeIPInp2d(lineN, background=optDefaults2d['IPbackground'], chifact=optDefaults2d['IPchifact'], useWeights=optDefaults2d['IPweights']):
    basedir = 'L%i'%lineN
    # Make inversion directories
    Inv = ['All', 'DipolePole', 'PoleDipole']
    # Inv = ['All']
    for inv in Inv:
        inpStr = ''
        inpStr += 'MESH FILE ../../mesh/mesh%i.txt\n'%lineN
        inpStr += 'COND FILE ../DC/dcinv2d.con\n'
        inpStr += 'OBS LOC_X ../../obs/L%i_IP_%s.obs\n'%(lineN, inv)
        inpStr += 'CHIFACT %.2f\n'%chifact
        inpStr += 'TOPO FILE ../../topo/topo%i.txt\n'%lineN
        inpStr += 'INIT_MOD VALUE %.2f\n'%background
        inpStr += 'REF_MOD VALUE %.2f\n'%background
        inpStr += 'ALPHA DEFAULT\n'
        inpStr += 'WAVE 2.5e-4 1.0 13\n'
        inpStr += 'STORE_ALL_MODELS TRUE\n'
        inpStr += 'INVMODE CG\n'
        inpStr += 'USE_MREF TRUE\n'
        inpStr += 'BOUNDS VALUE 0 100\n'
        if useWeights:
            inpStr += 'WEIGHT FILES ../../weight/S_weight.txt ../../weight/X_weight.txt ../../weight/Z_weight.txt\n'
        with open(basedir+'/'+inv+'/IP/ipinv2d.inp', 'w') as f: f.write(inpStr)
    
# Routines to read in UBC data

def readPredictedDatadc2d(fname,iterIsZero=False):
  if iterIsZero:
    df = pd.read_csv(fname,header=None,skiprows=1,usecols=[0,2,3,4,5],names=['ax','mx','nx','vpn','Appres'],delim_whitespace=True)
  else:
    df = pd.read_csv(fname,header=None,skiprows=1,usecols=[0,2,3,4],names=['ax','mx','nx','vpn'],delim_whitespace=True)
    
  df['ay'] = 0  
  df['my'] = 0
  df['ny'] = 0
  df['CP'] = 0.5*(df.mx+df.nx)
  df['Offset'] = np.abs(df.CP-df.ax)
  if not iterIsZero:
      df['Appres'] = appRes(df,'pdp')
      
  return df

def readObsDatadc2d(fname,hasErrs=True):
  if hasErrs:
    df = pd.read_csv(fname,header=None,skiprows=1,usecols=[0,2,3,4,5],names=['ax','mx','nx','vpn','Un'],delim_whitespace=True)
  else:
    df = pd.read_csv(fname,header=None,skiprows=1,usecols=[0,2,3,4],names=['ax','mx','nx','vpn'],delim_whitespace=True)
    
  df['ay'] = 0  
  df['my'] = 0
  df['ny'] = 0
  df['CP'] = 0.5*(df.mx+df.nx)
  df['Offset'] = np.abs(df.CP-df.ax)
  df['Appres'] = appRes(df,'pdp')
      
  return df
