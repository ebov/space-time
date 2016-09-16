import numpy as np
import pandas as pd
import matplotlib as mpl
from scipy import __version__ as scipy_version
from matplotlib.patches import Polygon

import re
import colorsys
import datetime as dt
import sys
import platform
import time
import math
from scipy.special import binom

def convertDate(x,start,end):
    """ Converts calendar dates between given formats """
    return dt.datetime.strftime(dt.datetime.strptime(x,start),end)

def decimalDate(date,fmt="%Y-%m-%d"):
    """ Converts calendar dates in specified format to decimal date. """
    adatetime=dt.datetime.strptime(date,fmt)
    year = adatetime.year
    boy = dt.datetime(year, 1, 1)
    eoy = dt.datetime(year + 1, 1, 1)
    return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds()))

#### code taken from Chris Slocum's website:
## http://schubert.atmos.colostate.edu/~cslocum/custom_cmap.html
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap

def hpd(data, level):
    """
    Return highest posterior density interval from a list,
    given the percent posterior density interval required.
    """
    d = list(data)
    d.sort()

    nData = len(data)
    nIn = int(round(level * nData))
    if nIn < 2 :
        return None
    #raise RuntimeError("Not enough data. N data: %s"%(len(data)))
 
    i = 0
    r = d[i+nIn-1] - d[i]
    for k in range(len(d) - (nIn - 1)) :
        rk = d[k+nIn-1] - d[k]
        if rk < r :
            r = rk
            i = k

    assert 0 <= i <= i+nIn-1 < len(d)
 
    return (d[i], d[i+nIn-1])

def overlap(a,b):
    """
    Return the elements shared by two lists in the following format:
    [overlap],[list 1 remainder],[list 2 remainder]
    """
    a_multiset = collections.Counter(a)
    b_multiset = collections.Counter(b)

    overlap = list((a_multiset & b_multiset).elements())
    a_remainder = list((a_multiset - b_multiset).elements())
    b_remainder = list((b_multiset - a_multiset).elements())

    return overlap, a_remainder, b_remainder

def column(data,col):
    return [row [col] for row in data]

def unique(o, idfun=repr):
    seen = {}
    return [seen.setdefault(idfun(e),e) for e in o if idfun(e) not in seen]

#### code stolen from seaborn to desaturate colours
def desaturate(color, prop):
    """Decrease the saturation channel of a color by some percent.
    Parameters
    ----------
    color : matplotlib color
        hex, rgb-tuple, or html color name
    prop : float
        saturation channel of color will be multiplied by this value
    Returns
    -------
    new_color : rgb tuple
        desaturated color code in RGB tuple representation
    """
    # Check inputs
    if not 0 <= prop <= 1:
        raise ValueError("prop must be between 0 and 1")

    # Get rgb tuple rep
    rgb = mpl.colors.colorConverter.to_rgb(color)

    # Convert to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)

    # Desaturate the saturation channel
    s *= prop

    # Convert back to rgb
    new_color = colorsys.hls_to_rgb(h, l, s)

    return new_color

def metricDistance(pointA,pointB):
    """ Calculate distance in kilometers along planet's surface from longitude and latitude. """
    R=6371 ## kilometers
    lon1,lat1=pointA
    lon2,lat2=pointB
    
    phi1=math.radians(lat1)
    phi2=math.radians(lat2)
    
    deltaPhi=math.radians(lat2-lat1)
    deltaLam=math.radians(lon2-lon1)
    
    a=math.sin(deltaPhi/2.0) * math.sin(deltaPhi/2.0) + math.cos(phi1) * math.cos(phi2) * math.sin(deltaLam/2.0) * math.sin(deltaLam/2.0)
    c = 2.0 * math.atan2(math.sqrt(a),math.sqrt(1-a))
    
    return R*c

## this function generates functions that are able to normalize custom value ranges to within a desired range
## it is used to normalize location coordinates within a country
def create_normalization(l,n_min,n_max):
    l_max=max(l)
    l_min=min(l)
    
    oldRange = (l_max-l_min)
    newRange = (n_max-n_min)
    
    def normalize(value):
        return (((value-l_min)*newRange)/float(oldRange)) + n_min
    return normalize

def Bernstein(n, k):
    """Bernstein polynomial.
    """
    coeff = binom(n, k)

    def _bpoly(x):
        return coeff * x ** k * (1 - x) ** (n - k)

    return _bpoly

def Bezier(points,start,end,num=10):
    """Build Bezier curve from points.
    """
    N = len(points)
    t = np.linspace(start, end, num)
    curve = np.zeros((num, 2))
    for ii in range(N):
        curve += np.outer(Bernstein(N - 1, ii)(t), points[ii])
        
    return curve

def Bezier_control(pointA,pointB,height):
    """ 
    Given a line defined by 2 points A & B, 
    find a third point at a given distance that defines a line perpendicular to line AB which intercepts AB at its midpoint.
    Equation derived by Luiz Max Fagundes de Carvalho (University of Edinburgh).
    """
    x1,y1=pointA
    x2,y2=pointB
    
    sign=1
    if x1>x2:
        sign=-1

    slope = (y2-y1) / (x2-x1)
    d=np.sqrt((y2-y1)**2 + (x2-x1)**2)
        
    h=np.sqrt(height**2+(d**2)/4.0)
    
    n1=x1+h*np.cos(np.arctan(2*height/float(d))+np.arctan(slope))*sign
    n2=y1+h*np.sin(np.arctan(2*height/float(d))+np.arctan(slope))*sign

    return (n1,n2)

## path to master directory
path_to_dropbox='/Users/evogytis/Dropbox/Ebolavirus_Phylogeography/'

## import location information
map_path=path_to_dropbox+'Data_Processing/June2016/'
filenames=['location_data_v3.geojson']

## path to population centroids
loc_path=path_to_dropbox+'Sequences/locations.txt'

## normalized location coordinates after PCA
PCA_path=path_to_dropbox+'plotting_assists/EBOV_location_PCA.txt'

## path to case numebrs
cases_path=path_to_dropbox+'Case counts by district/EBOV_maxCases.csv'

## define a path where figures will be saved
local_output='/Users/evogytis/Documents/EBOV_output/'

## contains the coordinates of the international borders
border_path=path_to_dropbox+'plotting_assists/international_border.txt'

## actual names of locations (with diacritics)
standard=path_to_dropbox+'Maps/standardDistricts.tsv'

## tinkering with details of where location names should be plotted
textpath=path_to_dropbox+'plotting_assists/EBOV_plotting_assist.txt'

## path to where location border sharing information is
border_sharing=path_to_dropbox+'Data_Processing/June2016/shared_international_border_81.txt'

## set location colouring based on PCA1 coordinate
normalize_by='PCA1'

global xlimits
xlimits=[]

global ylimits
ylimits=[]

global required_countries
required_countries=[]

def setFocusCountries(focus=['SLE','LBR','GIN','SEN','GNB','CIV','MLI']):
	
	for country in focus:
		required_countries.append(country)
	
	if len(focus)>3:
		xlimits.append(-15.5)
		xlimits.append(-6.55)
		ylimits.append(4.3)
		ylimits.append(14.5)
	else:
		xlimits.append(-15.5)
		xlimits.append(-7.3)
		ylimits.append(4.3)
		ylimits.append(12.7)

def status():
	print '%-20s%10s (%s)'%('Operating system:',platform.system(),platform.release())
	print '%-20s%10s'%('Python version:','.'.join(map(str,sys.version_info[:3])))
	print '%-20s%10s'%('Numpy version:',np.__version__)
	print '%-20s%10s'%('Pandas version:',pd.__version__)
	print '%-20s%10s'%('matplotlib version:',mpl.__version__)
	print '%-20s%10s'%('scipy version:',scipy_version)
	print '\nThis notebook was last run on:\n%s\t%s'%(dt.datetime.strftime(dt.datetime.now(),'%A\t%Y-%b-%d\t%H:%M'),time.tzname[time.daylight])

global colours ## dict with country colours
colours={}

def setColourMaps():
	## Start loading data and creating various plotting assists
	## create new colour maps for each country with desaturated colours
	desaturation=0.8
	initial_colours={'SLE':[(234,238,246),(141,160,203),(57,77,124),(33,44,67)],
					 'LBR':[(254,233,226),(252,141,98),(177,52,4),(97,31,3)],
					 'GIN':[(234,246,242),(102,194,165),(52,129,104),(31,69,57)],
					 'SEN':[(240,240,240),(175,175,175),(90,90,90),(50,50,50)],
					 'MLI':[(240,240,240),(175,175,175),(90,90,90),(50,50,50)],
					 'GNB':[(240,240,240),(175,175,175),(90,90,90),(50,50,50)],
					 'CIV':[(240,240,240),(175,175,175),(90,90,90),(50,50,50)],
					 '?':[(240,240,240),(175,175,175),(90,90,90),(50,50,50)]}

	for country in initial_colours.keys():
		colours[country]=make_cmap([desaturate(np.array(x)/256.0,desaturation) for x in initial_colours[country]])

global locations ## list of administrative region names
locations=[]

global polygons ## dict of polygons
polygons={}

global location_points ## dict of polygon coordinates
location_points={}

global location_to_country ## maps location to its country
location_to_country={}

global map_to_actual ## dict with actual accented names of prefectures
map_to_actual={}

global textCorrection ## dict indicating how location labels should be aligned
textCorrection={}

global cases_byCountry ## dict of case numbers over time by country
cases_byCountry={}

global cases_byLocation ## dict of case numbers over time by location
cases_byLocation={}

global PCA ## dict of PCA coordinates of each location
PCA={}

global popCentres ## dict of population centroid coordinates
popCentres={}

global normalized_coords ## dict of normalized location coordinates
normalized_coords={}

global dates ## list of epiweeks
dates=[]

global maxByCountry ## dict of country to cumulative cases
maxByCountry={country:0 for country in required_countries}

global global_border ## list of international border fragments
global_border=[]

global shared_border ## boolean dictionary of border sharing between locations
shared_border={}

def loadData():

	for fname in filenames: ## iterate through country files
		handle=pd.json.load(open(map_path+fname,'r')) ## load data
		features=handle['features']

		for loc in features: ## iterate through features (locations)
			poly = np.array(loc['geometry']['coordinates']) ## get coordinates
			country=loc['properties']['ISO']
			location=loc['properties']['location'] ## standardised location name
		
			location_to_country[location]=country
			locations.append(location)
		
			polygons[location]=[]
			location_points[location]=[]
			if loc['geometry']['type']=='MultiPolygon': ## multiple parts detected
				for part in poly:
					for coords in part:
						xs=column(coords,0)
						ys=column(coords,1)
						location_points[location].append(np.vstack(zip(xs,ys)))
			if loc['geometry']['type']=='Polygon': ## location is single part
				for coords in poly:
					xs=column(coords,0)
					ys=column(coords,1)
					location_points[location].append(np.vstack(zip(xs,ys)))
		
			complete_location=[]
			for part in location_points[location]:
				complete_location.append(Polygon(part,True))

			polygons[location]=complete_location
	location_to_country['WesternArea']='SLE'
	
	
	## Load actual names of locations
	for line in open(standard,'r'):
		l=line.strip('\n').split('\t')
		## map parts of a district name to its standard name
		# split by space, keep last word in lower case
		# e.g River Gee = gee, Grand Cape Mount = mount.
		actual=l[1]
		if 'actual' not in actual:
			actual=actual.decode('utf8')
			map_to_actual[l[-1]]=actual

	## Load text tinkering
	for line in open(textpath,'r'):
		location,t,r=line.strip('\n').split('\t')
		textCorrection[location]=(int(t),int(r))

	## Load case numbers
	for line in open(cases_path,'r'):
		l=line.strip('\n').split(',')
		if l[0]=='country':
			for d in l[3:]:
				dates.append(convertDate(d,'%Y-%b-%d','%Y-%m-%d'))
		else:
			if cases_byCountry.has_key('%s'%(l[0]))==False:
				cases_byCountry['%s'%(l[0])]={}
			
			for i,x in enumerate(dates):
				if cases_byCountry['%s'%(l[0])].has_key(x)==False:
					cases_byCountry['%s'%(l[0])][x]=0
				
				if l[3+i]!='':
					cases_byCountry['%s'%(l[0])][x]+=int(l[3+i])
				else:
					cases_byCountry['%s'%(l[0])][x]+=0
				
				
			d='%s'%(l[2])
			if cases_byLocation.has_key('%s'%(d))==False:
				cases_byLocation['%s'%(d)]={}
			
			for i,x in enumerate(dates):
				if cases_byLocation['%s'%(d)].has_key(x)==False:
					cases_byLocation['%s'%(d)][x]=0
				
				if l[3+i]!='':
					cases_byLocation['%s'%(d)][x]+=int(l[3+i])
				else:
					cases_byLocation['%s'%(d)][x]+=0
				
	for loc in locations:
		if cases_byLocation.has_key(loc)==False:
			cases_byLocation[loc]={d:0 for d in dates}
			   
	totalCaseCounts={x:sum(cases_byCountry[x].values()) for x in cases_byCountry.keys()}

	# for location in cases_byLocation.keys():
# 		country=location_to_country[location]
# 		cases_in_location=sum(cases_byLocation[location].values())
# 		if maxByCountry[country]<=cases_in_district:
# 			maxByCountry[country]=cases_in_district

	## Load PCA coordinates of population centroids
	for line in open(PCA_path,'r'):
		l=line.strip('\n').split('\t')
		if l[0]=='location':
			pass
		else:
			PCA[l[0]]=(float(l[-2]),float(l[-1]))

	## PCA coordinates already normalized
	normalization=create_normalization([0,1],0.2,1.0)

	for loc in locations:
		if normalize_by=='PCA1':
			normalized_coords[loc]=normalization(PCA[loc][0])
		elif normalize_by=='PCA2':
			normalized_coords[loc]=normalization(PCA[loc][1])

	fix={'WesternArea':('WesternRural','WesternUrban')}

	for f in fix.keys():
		being_fixed=f
		fixed_with=fix[being_fixed]
	
		aggregate=[]
	
		for f2 in fixed_with:
			if normalize_by=='PCA1':
				aggregate.append(PCA[f2][0])
			elif normalize_by=='PCA2':
				aggregate.append(PCA[f2][1])
			
		normalized_coords[being_fixed]=normalization(np.mean(aggregate))
	
	## import population centroids
	for line in open(loc_path,'r'):
		if 'Location' not in line:
			location,lon,lat=line.strip('\n').split('\t')
			popCentres[location]=(float(lon),float(lat))
			
	## import coordinates of the international borders
	for line in open(border_path,'r'):
		cerberus=re.findall('\([0-9\.\- ,]+\)',line.strip('\n'))
		global_border.append([map(float,x[1:-1].split(',')) for x in cerberus])
		
	## set up dictionary of border sharing
	for line in open(border_sharing,'r'):
	    l=line.strip('\n').split('\t')
	    if l[0]=='locA':
	        pass
	    else:
	        locA,locB,share,international=l
	        if int(share)==0:
	            s=False
	        else:
	            s=True
	        
	        if shared_border.has_key(locA)==False:
	            shared_border[locA]={locB:s}
	            
	        if shared_border.has_key(locB)==False:
	            shared_border[locB]={locA:s}
	            
	        shared_border[locA][locB]=s
	        shared_border[locB][locA]=s
	    