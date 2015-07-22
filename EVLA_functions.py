######################################################################
# Functions for EVLA pipeline script
# Runs on CASA 3.3.0
# 02/14/12 CJC
######################################################################

# From ALMA/Nick Elias:

# Usage:
# out = selectReferenceAntenna( '<ms>', '<spws>', '<fields>' )

# ------------------------------------------------------------------------------

# selectReferenceAntenna

# Description:
# ------------
# This function selects the reference antenna based on x-y position within the
# the array as well as the number of completely flagged rows.

# NB: This function is meant to be a temporary stop gap.  The pipeline version
# will be organized differently and perhaps include additional criteria.

# NB: This function determines the reference antenna by separate distance and
# flag "scores" and creating from them a combined "score."

# NB: This function does not use intents to determine which fields are
# calibrators.  The user must select fields.  

# NB: This function does look for phase jumps.

# Inputs:
# -------
# ms     - The MS name, specified as a python string.  If the MS is not in the
#          present working directory, the relative or absolute paths must be
#          specified.
# spws   - The spectral window IDs, specified as a comma-delimited python string.
#          The default is '' (all spectral window IDs).  NB: Do not use names
#          from the NAME column in the SPECTRAL_WINDOW subtable, use only the
#          numerical IDs.  Also, do no include channel selection.
# fields - The field IDs, specified as a comma-delimited python string.  The
#          default is '' (all field IDs).  NB: Do not use names from the NAME
#          column in the FIELD subtable.  Use only the numerical IDs.

# Outputs:
# --------
# A python combination, returned via the function value, whose elements are:
# out[0] - A python string containing the name of the reference antenna ID (not
#          the numerical antenna IDs).
# out[1] - A numpy array of strings containing the antenna IDs, ordered
#          according to total score (the first element is the reference antenna
#          ID).
# out[2] - A python dictionary containing the total scores.  The keywords are the
#          antenna IDs and the elements are the total scores.
# out[3] - A python dictionary containing the distance scores.  The keywords are
#          the antenna IDs and the elements are the distance scores.
# out[4] - A python dictionary containing the flag row scores.  The keywords are
#          the antenna IDs and the elements are the flag row scores.

# Modification history:
# ---------------------
# 2012 Feb 13 - Nick Elias, NRAO
#               This function was pilfered from the original pipeline code and
#               cleaned up.  It was converted from a member function to a global
#               function.  All references to member functions were removed.  The
#               spectral window and field inputs are now comma-delimited strings
#               with '' defaults (before, they were lists of integers with no
#               defaults).

# ------------------------------------------------------------------------------

def selectReferenceAntenna( ms, spws='', fields='' ):

	# Create the local instance of the table, measures, and quanta tools

	tbLoc = casa.__tablehome__.create()
	meLoc = casa.__measureshome__.create()
	qaLoc = casa.__quantahome__.create()


	# Fix the spectral window input

	if spws != '':

		spwSList = spws.split( ',' )

		spwList = []
		for s in spwSList: spwList.append( int(s) )

	else:

		tbLoc.open( ms + '/SPECTRAL_WINDOW' )
		rSPW = range( tbLoc.nrows() )
		tbLoc.close()

		spwList = []
		for s in rSPW: spwList.append( int(s) )


	# Fix the field input

	if fields != '':

		fieldSList = fields.split( ',' )

		fieldList = []
		for f in fieldSList: fieldList.append( int(f) )

	else:

		tbLoc.open( ms + '/FIELD' )
		rField = range( tbLoc.nrows() )
		tbLoc.close()

		fieldList = []
		for f in rField: fieldList.append( int(f) )


	# Read the antenna positions.  These seem to be in metres rel to centre
	# of the Earth.
     
	tbLoc.open( '%s/ANTENNA' % ms )

	position = tbLoc.getcol( 'POSITION' )
	flag_row = tbLoc.getcol( 'FLAG_ROW' )
	name = tbLoc.getcol( 'NAME' )
	position_keywords = tbLoc.getcolkeywords( 'POSITION' )

	tbLoc.close()


	# Make a position Measure for each antenna, this stores info in
	# terms of long, lat and distance from centre of Earth.

	antenna_position = {}

	for row,ant in enumerate( name ):
		if not flag_row[row]:
			rf = position_keywords['MEASINFO']['Ref']
			v0 = qaLoc.quantity( position[0,row],
			    position_keywords['QuantumUnits'][0] )
			v1 = qaLoc.quantity( position[1,row],
			    position_keywords['QuantumUnits'][1] )
			v2 = qaLoc.quantity( position[2,row],
			    position_keywords['QuantumUnits'][2] )
			ant_pos = meLoc.position( rf=rf, v0=v0, v1=v1, v2=v2 )
			antenna_position[ant] = ant_pos


	# Store the longs and lats of the antennas in lists - convert to
	# canonical units.

	longs = {}
	lats = {}
	radii = {}

	for ant in name:

		position = antenna_position[ant]
		radius = position['m2']['value']
		radius_unit = position['m2']['unit']
		radius_quantum = qaLoc.quantity( radius, radius_unit )
		radius_quantum = qaLoc.convert( radius_quantum, 'm' )
		radius = qaLoc.getvalue( radius_quantum )

		long = position['m0']['value']
		long_unit = position['m0']['unit']
		long_quantum = qaLoc.quantity( long, long_unit )
		long_quantum = qaLoc.convert( long_quantum, 'rad' )
		long = qaLoc.getvalue( long_quantum )

		lat = position['m1']['value']
		lat_unit = position['m1']['unit']
		lat_quantum = qaLoc.quantity( lat, lat_unit )
		lat_quantum = qaLoc.convert( lat_quantum, 'rad' )
		lat = qaLoc.getvalue( lat_quantum )

		longs[ant] = long
		lats[ant] = lat
		radii[ant] = radius


	# Get a median 'centre' for the array and derive x, y relative
	# to that.

	long_vals = np.array( longs.values() )
	long_vals -= np.median( long_vals )

	lat_vals = np.array( lats.values() )
	radii_vals = np.array( radii.values() )

	names = longs.keys()


	# Multiply longs by cos(lat) and radius to convert to metres.

	x = long_vals * np.cos(lat_vals) * radii_vals
	y = lat_vals * radii_vals


	# Make x,y relative to 'centre' of array.

	x -= np.median(x)
	y -= np.median(y)

	antenna_distance = {}

	for i,ant in enumerate(names):
		antenna_distance[ant] = np.sqrt( pow(x[i],2) + pow(y[i],2) )

	antennaInfoRead = True


        # Assign points to antennas based on their distance from the array
        # centre.

	distance_vals = antenna_distance.values()
	names = antenna_distance.keys()
	distance_argsort = np.argsort( np.array(distance_vals) )
	distance_points = {}

	for k,i in enumerate( distance_argsort ):
		ant = names[i]
		distance_points[ant] = len(distance_argsort) - k


	# Get the number of unflagged rows for each antenna in these spectral
	# windows / fields.

	tbLoc.open( '%s' % ms )

        taql = '(FIELD_ID IN %s)' % fieldList
	taql += ' && (DATA_DESC_ID IN %s)' % spwList
	taql += ' && (NOT(FLAG_ROW))'

	nrows = {}
	ids = np.arange( len(name) )

	for ant in name:
		ant_id = ids[name==ant][0]
		subTable = tbLoc.query(
		  '%s && (ANTENNA1==%s || ANTENNA2==%s)' % (taql,ant_id,ant_id)
		)
		nrows[ant] = subTable.nrows()

        tbLoc.close()


        # Assign points to antennas for the amount of data they have.

        nrows_vals = np.array( nrows.values(), np.float )
        names = np.array( nrows.keys() )
        nrows_normalised = (nrows_vals / max(nrows_vals)) * len(nrows_vals)
        nrows_points_array = map( int, nrows_normalised )

        nrows_points = {}

        for i,ant in enumerate(names):
		nrows_points[ant] = nrows_points_array[i]


        # Add up points for antennas.

	total_points = {}

	for ant in name:
		total_points[ant] = distance_points[ant] + nrows_points[ant]


        # Select antenna with the highest total score.

	total_points_array = np.array( total_points.values() )
	total_points_array = - total_points_array
	merit_argsort = np.argsort( total_points_array )

	sorted_antennas = names[merit_argsort]
	reference_antenna = names[merit_argsort[0]]


	# Return the reference antenna, sorted antenna list, and points

	return reference_antenna, sorted_antennas, total_points, \
	    distance_points, nrows_points

######################################################################

def uniq(inlist):
   uniques = []
   for item in inlist:
      if item not in uniques:
         uniques.append(item)
   return uniques

######################################################################

def find_EVLA_band(frequency):
   FLOW = [ 0.0e6, 150.0e6, 700.0e6, 2.0e9, 4.0e9, 8.0e9, 12.0e9, 18.0e9, 26.5e9, 40.0e9 ]
   FHIGH = [ 150.0e6, 700.0e6, 2.0e9, 4.0e9, 8.0e9, 12.0e9, 18.0e9, 26.5e9, 40.0e9, 56.0e9 ]
   BBAND = [ '4', 'P', 'L', 'S', 'C', 'X', 'U', 'K', 'A', 'Q' ]

   band = '?'
   for ii in range(0,len(FLOW)):
      if ((frequency > FLOW[ii]) and (frequency < FHIGH[ii])):
         band = BBAND[ii]
   return band

######################################################################

def find_standards(positions):
# set the max separation as ~1'
   MAX_SEPARATION = 60*2.0e-5
   position_3C48 = me.direction('j2000', '1h37m41.299', '33d9m35.133')
   fields_3C48 = []
   position_3C138 = me.direction('j2000', '5h21m9.886', '16d38m22.051')
   fields_3C138 = []
   position_3C147 = me.direction('j2000', '5h42m36.138', '49d51m7.234')
   fields_3C147 = []
   position_3C286 = me.direction('j2000', '13h31m8.288', '30d30m23.959')
   fields_3C286 = []

   for ii in range(0,len(positions)):
      position = me.direction('j2000', str(positions[ii][0])+'rad', str(positions[ii][1])+'rad')
      separation = me.separation(position,position_3C48)['value'] * pi/180.0
      if (separation < MAX_SEPARATION):
         fields_3C48.append(ii)
      else:
         separation = me.separation(position,position_3C138)['value'] * pi/180.0
         if (separation < MAX_SEPARATION):
            fields_3C138.append(ii)
         else:
            separation = me.separation(position,position_3C147)['value'] * pi/180.0
            if (separation < MAX_SEPARATION):
               fields_3C147.append(ii)
            else:
               separation = me.separation(position,position_3C286)['value'] * pi/180.0
               if (separation < MAX_SEPARATION):
                  fields_3C286.append(ii)

   fields = [ fields_3C48, fields_3C138, fields_3C147, fields_3C286 ]
   return fields

######################################################################

def getFlaggedSoln(calTable):
    '''
    This method will look at the specified calibration table and return the
    fraction of flagged solutions for each Antenna, SPW, Poln.  This assumes
    that the specified cal table will not have any channel dependent flagging.

    return structure is a dictionary with AntennaID and Spectral Window ID
    as the keys and returns a list of fractional flagging per polarization in
    the order specified in the Cal Table.

    Example:
    execfile('getFlaggedSoln.py')
    result = getFlaggedSoln('testdelay.k')
    result[17][4] = [1.0, 1.0]
    '''

    from taskinit import tbtool
    mytb = tbtool.create()

    mytb.open(calTable)
    antCol = mytb.getcol('ANTENNA1')
    flagCol = mytb.getcol('FLAG')
    calDescCol = mytb.getcol('CAL_DESC_ID')
    mytb.close()

    mytb.open(calTable+"/CAL_DESC")
    spwCol = mytb.getcol('SPECTRAL_WINDOW_ID')[0]
    mytb.done()
    
    # Initialize a list to hold the results
    # The order here is row[antenna][calDescCol]([polarization] for flags)
    rows = [[0] * (max(calDescCol)+1) for idx in range(max(antCol)+1)]
    flags = [[[0] * len(flagCol) for idx in range(max(calDescCol)+1)]
             for idx2 in range(max(antCol)+1)]

    # Ok now go through and for each antenna and calDesc get the number
    # of flagged and total number of entries
    # Here we assume that there is no channel dependent information
    for idx in range(len(antCol)):
        rows[antCol[idx]][calDescCol[idx]] += 1
        for poln in range(len(flagCol)):
            if flagCol[poln][0][idx]:
                flags[antCol[idx]][calDescCol[idx]][poln] += 1
                
    # Now create the output dictionary
    outDict = {}
    for antIdx in range(max(antCol)+1):
        outDict[antIdx] = {}
        for calDescIdx in range(max(calDescCol)+1):
            if rows[antIdx][calDescIdx] > 0:
                outDict[antIdx][spwCol[calDescIdx]] = \
                       [pF/float(rows[antIdx][calDescIdx]) \
                        for pF in flags[antIdx][calDescIdx]]



    return outDict
    
######################################################################

import math
def getOptimumSize(size):
    '''
    This method takes as input the a size parameter.  The return is the smallest
    integer Y which satisfies the following conditions:
    * Y > size
    * Y = 2^a*3^b*5^c where a,b, and c are non-negative integers and at least one
    of a or b is 0 and c is nonzero
    '''
    def evaluate(pow2, pow3, pow5):
        # Convience method to calculate the value given multiples
        return int(math.pow(2,pow2) *math.pow(3,pow3)*math.pow(5,pow5))
    
    max5 = int(math.ceil(math.log(size,5)))
    returnValue = evaluate(0, 0, max5)
    for pow5 in range(max5,0,-1):
        pow2 = math.ceil(math.log(size/math.pow(5,pow5),2))
        if not pow2 < 0:
            returnValue = min(returnValue, evaluate(pow2,0,pow5))

        pow3 = math.ceil(math.log(size/math.pow(5,pow5),3))
        if not pow3 < 0:
            returnValue = min(returnValue, evaluate(0,pow3,pow5))
    return returnValue

######################################################################

import urllib2
import datetime

def correct_ant_posns (vis_name, print_offsets=False):
    '''
    Given an input visibility MS name (vis_name), find the antenna
    position offsets that should be applied.  This application should
    be via the gencal task, using caltype='antpos'.

    If the print_offsets parameter is True, will print out each of
    the found offsets (or indicate that none were found), otherwise
    runs silently.

    A list is returned where the first element is the returned error
    code, the second element is a string of the antennas, and the 
    third element is a list of antenna Bx,By,Bz offsets.  An example 
    return list might look like:
    [ 0, 'ea01,ea19', [0.0184, -0.0065, 0.005, 0.0365, -0.0435, 0.0543] ]

    Usage examples:

       CASA <1>: antenna_offsets = correct_ant_posns('test.ms')
       CASA <2>: if (antenna_offsets[0] == 0):
       CASA <3>:     gencal(vis='test.ms', caltable='cal.G', \
                     caltype='antpos', antenna=antenna_offsets[1], \
                     parameter=antenna_offsets[2])

    This function does NOT work for VLA datasets, only EVLA.  If an
    attempt is made to use the function for VLA data (prior to 2010),
    an error code of 1 is returned.

    The offsets are retrieved over the internet.  A description and the
    ability to manually examine and retrieve offsets is at:
    http://www.vla.nrao.edu/astro/archive/baselines/
    If the attempt to establish the internet connection fails, an error
    code of 2 is returned.

    Uses the same algorithm that the AIPS task VLANT does.


    bjb
    nrao
    spring 2012
    '''

    MONTHS = [ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' ]
    URL_BASE = 'http://www.vla.nrao.edu/cgi-bin/evlais_blines.cgi?Year='

#
# get start date+time of observation
#
    observation = tb.open(vis_name+'/OBSERVATION')
    time_range = tb.getcol('TIME_RANGE')
    tb.close()
    MJD_start_time = time_range[0][0] / 86400
    q1 = qa.quantity(time_range[0][0],'s')
    date_time = qa.time(q1,form='ymd')
# date_time looks like: '2011/08/10/06:56:49'
    [obs_year,obs_month,obs_day,obs_time_string] = date_time.split('/')
    if (int(obs_year) < 2010):
        if (print_offsets):
            print 'Does not work for VLA observations'
        return [1, '', []]
    [obs_hour,obs_minute,obs_second] = obs_time_string.split(':')
    obs_time = 10000*int(obs_year) + 100*int(obs_month) + int(obs_day) + \
               int(obs_hour)/24.0 + int(obs_minute)/1440.0 + \
               int(obs_second)/86400.0

#
# get antenna to station mappings
#
    observation = tb.open(vis_name+'/ANTENNA')
    ant_names = tb.getcol('NAME')
    ant_stations = tb.getcol('STATION')
    tb.close()
    ant_num_stas = []
    for ii in range(len(ant_names)):
        ant_num_stas.append([int(ant_names[ii][2:]), ant_names[ii], \
                            ant_stations[ii], 0.0, 0.0, 0.0, False])

    correction_lines = []
    current_year = datetime.datetime.now().year
# first, see if the internet connection is possible
    try:
        response = urllib2.urlopen(URL_BASE + '2010')
    except URLError, err:
        if (print_offsets):
            print 'No internet connection to antenna position correction URL ', \
                  err.reason
        return [2, '', []]
    response.close()
    for year in range(2010,current_year+1):
        response = urllib2.urlopen(URL_BASE + str(year))
        html = response.read()
        response.close()
        html_lines = html.split('\n')
        for correction_line in html_lines:
            for month in MONTHS:
                if (correction_line.find(month) >= 0):
                    correction_lines.append(str(year)+' '+correction_line)
                    break

    corrections_list = []
    for correction_line in correction_lines:
        correction_line_fields = correction_line.split()
        if (len(correction_line_fields) > 9):
            [c_year, moved_date, obs_date, put_date, put_time_str, ant, pad, Bx, By, Bz] = correction_line_fields
            s_moved = moved_date[:3]
            i_month = 1
            for month in MONTHS:
                if (moved_date.find(month) >= 0):
                    break
                i_month = i_month + 1
            moved_time = 10000 * int(c_year) + 100 * i_month + \
                         int(moved_date[3:])
        else:
            [c_year, obs_date, put_date, put_time_str, ant, pad, Bx, By, Bz] = correction_line_fields
            moved_date = '     '
            moved_time = 0
        s_obs = obs_date[:3]
        i_month = 1
        for month in MONTHS:
            if (s_obs.find(month) >= 0):
                break
            i_month = i_month + 1
        obs_time_2 = 10000 * int(c_year) + 100 * i_month + int(obs_date[3:])
        s_put = put_date[:3]
        i_month = 1
        for month in MONTHS:
            if (s_put.find(month) >= 0):
                break
            i_month = i_month + 1
        put_time = 10000 * int(c_year) + 100 * i_month + int(put_date[3:])
        [put_hr, put_min] = put_time_str.split(':')
        put_time += (int(put_hr)/24.0 + int(put_min)/1440.0)
        corrections_list.append([c_year, moved_date, moved_time, obs_date, obs_time_2, put_date, put_time, int(ant), pad, float(Bx), float(By), float(Bz)])

    for correction_list in corrections_list:
        [c_year, moved_date, moved_time, obs_date, obs_time_2, put_date, put_time, ant, pad, Bx, By, Bz] = correction_list
        ant_ind = -1
        for ii in range(len(ant_num_stas)):
            ant_num_sta = ant_num_stas[ii]
            if (ant == ant_num_sta[0]):
                ant_ind = ii
                break
        if ((ant_ind == -1) or (ant_num_sta[6])):
# the antenna in this correction isn't in the observation, or is done, 
# so skip it
            pass
        ant_num_sta = ant_num_stas[ant_ind]
        if (moved_time):
# the antenna moved
            if (moved_time > obs_time):
# we are done considering this antenna
                ant_num_sta[6] = True
            else:
# otherwise, it moved, so the offsets should be reset
                ant_num_sta[3] = 0.0
                ant_num_sta[4] = 0.0
                ant_num_sta[5] = 0.0
        if ((put_time > obs_time) and (not ant_num_sta[6]) and \
            (pad == ant_num_sta[2])):
# it's the right antenna/pad; add the offsets to those already accumulated
            ant_num_sta[3] += Bx
            ant_num_sta[4] += By
            ant_num_sta[5] += Bz

    ants = []
    parms = []
    for ii in range(len(ant_num_stas)):
        ant_num_sta = ant_num_stas[ii]
        if ((ant_num_sta[3] != 0.0) or (ant_num_sta[4] != 0.0) or \
            (ant_num_sta[3] != 0.0)):
            if (print_offsets):
                print "offsets for antenna %4s : %8.5f  %8.5f  %8.5f" % \
                      (ant_num_sta[1], ant_num_sta[3], ant_num_sta[4], ant_num_sta[5])
            ants.append(ant_num_sta[1])
            parms.append(ant_num_sta[3])
            parms.append(ant_num_sta[4])
            parms.append(ant_num_sta[5])
    if ((len(parms) == 0) and print_offsets):
        print "No offsets found for this MS"
    ant_string = ','.join(["%s" % ii for ii in ants])
    return [ 0, ant_string, parms ]

def logprint(msg):
    print (msg)
    outlog.write(msg+"\n")
    outlog.flush()
