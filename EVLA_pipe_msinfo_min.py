######################################################################

# GET SOME INFORMATION FROM THE MS THAT WILL BE NEEDED LATER

outlog=open(fluxlogname,'w')
outlog.close()

# Identify spw information

tb.open(msname+'/SPECTRAL_WINDOW')
channels = tb.getcol('NUM_CHAN')
originalBBClist = tb.getcol('BBC_NO')
spw_bandwidths=tb.getcol('TOTAL_BANDWIDTH')
reference_frequencies = tb.getcol('REF_FREQUENCY')
center_frequencies = []
for ii in range(len(reference_frequencies)):
    center_frequencies.append(reference_frequencies[ii]+spw_bandwidths[ii]/2)
tb.close()

numSpws = len(channels)

# Set up spw selection for initial gain solutions

tst_delay_spw=''

for ispw in range(numSpws):
   endch1=int(channels[ispw]/3.0)
   endch2=int(2.0*channels[ispw]/3.0)+1
   if (ispw<max(range(numSpws))):
      tst_delay_spw=tst_delay_spw+str(ispw)+':'+str(endch1)+'~'+str(endch2)+','
   else:
      tst_delay_spw=tst_delay_spw+str(ispw)+':'+str(endch1)+'~'+str(endch2)

tst_bpass_spw=tst_delay_spw

# Identify number of fields, positions, and source IDs

tb.open(msname+'/FIELD')
sources = tb.getcol('NAME')
numFields = tb.nrows()
field_positions = tb.getcol('PHASE_DIR')
source_ids = tb.getcol('SOURCE_ID')
tb.close()

# Map source names to spws

tb.open(msname+'/SOURCE')
num_sources = tb.nrows()
source_spws = []
for ii in range(0,num_sources):
   subtable = tb.query('SOURCE_ID=%s'%ii)
   source_spws.append(list(np.unique(subtable.getcol('SPECTRAL_WINDOW_ID'))))
tb.close()

# Identify scan numbers, and map scans to field ID

tb.open(msname)
scanNums = sorted(np.unique(tb.getcol('SCAN_NUMBER')))
tb.close()

tb.open(msname)
field_scans = []
for ii in range(0,numFields):
   subtable = tb.query('FIELD_ID==%s'%ii)
   field_scans.append(list(np.unique(subtable.getcol('SCAN_NUMBER'))))
tb.close()

## field_scans is now a list of lists containing the scans for each field.
## so, to access all the scans for the fields, you'd:
#
#for ii in range(0,len(field_scans)):
#   for jj in range(0,len(field_scans[ii]))
#
## the jj'th scan of the ii'th field is in field_scans[ii][jj]

# Identify intents

tb.open(msname+'/STATE')
intents=tb.getcol('OBS_MODE')
tb.close()

# Figure out integration time used

ms.open(msname)
scan_summary = ms.getscansummary()
ms_summary = ms.summary()
ms.close()
startdate=float(ms_summary['header']['BeginTime'])
#
# scan list
#
integ_scan_list = []
for scan in scan_summary['summary']:
    integ_scan_list.append(int(scan))
sorted_scan_list = sorted(integ_scan_list)
#
# find max and median integration times
#
integration_times = []
spws_per_scan = {}
for ii in sorted_scan_list:
   integration_times.append(scan_summary['summary'][str(ii)]['0']['IntegrationTime'])
   spws_per_scan[str(ii)] = list(scan_summary['summary'][str(ii)]['0']['SpwIds'])

maximum_integration_time = max(integration_times)
median_integration_time = np.median(integration_times)

int_time=maximum_integration_time

# Find scans for quacking

scan_list = [1]
old_scan = scan_summary['summary'][str(sorted_scan_list[0])]['0']
old_field = old_scan['FieldId']
old_spws = old_scan['SpwIds']
for ii in range(1,len(sorted_scan_list)):
    new_scan = scan_summary['summary'][str(sorted_scan_list[ii])]['0']
    new_field = new_scan['FieldId']
    new_spws = new_scan['SpwIds']
    if ((new_field != old_field) or (set(new_spws) != set(old_spws))):
        scan_list.append(sorted_scan_list[ii])
        old_field = new_field
        old_spws = new_spws
quack_scan_string = ','.join(["%s" % ii for ii in scan_list])

# For 1 GHz wide basebands, figure out which spws are associated
# with the edges of the baseband filters

sorted_frequencies = sorted(reference_frequencies)
sorted_indices = []

for ii in range (0,len(sorted_frequencies)):
    for jj in range (0,len(reference_frequencies)):
        if (sorted_frequencies[ii] == reference_frequencies[jj]):
            sorted_indices.append(jj)

spwList = []
BBC_bandwidths = []
ii = 0

while (ii < len(sorted_frequencies)):
    upper_frequency = sorted_frequencies[ii] + spw_bandwidths[sorted_indices[ii]]
    BBC_bandwidth = spw_bandwidths[sorted_indices[ii]]
    thisSpwList = [sorted_indices[ii]]
    jj = ii + 1
    while (jj < len(sorted_frequencies)):
        lower_frequency = sorted_frequencies[jj]
        if ((fabs(lower_frequency - upper_frequency) < 1.0) and \
            (originalBBClist[sorted_indices[ii]] == originalBBClist[sorted_indices[jj]])):
            thisSpwList.append(sorted_indices[jj])
            upper_frequency += spw_bandwidths[sorted_indices[jj]]
            BBC_bandwidth += spw_bandwidths[sorted_indices[jj]]
            jj += 1
            ii += 1
        else:
            jj = len(sorted_frequencies)
    spwList.append(thisSpwList)
    BBC_bandwidths.append(BBC_bandwidth)
    ii += 1

# spwList is now a list of lists of contiguous spws, which have
# bandwidths in BBC_bandwidths

low_spws = []
high_spws = []

for ii in range(0,len(BBC_bandwidths)):
   if (BBC_bandwidths[ii] > 1.0e9):
      low_spws.append(spwList[ii][0])
      high_spws.append(spwList[ii][len(spwList[ii])-1])

if (len(low_spws) != len(high_spws)):
   raise Exception("Error! Something is wrong with the spw identification")

# Identify scans and fields associated with different calibrator intents

if 1:
   bandpass_state_IDs = []
   delay_state_IDs = []
   flux_state_IDs = []
   polarization_state_IDs = []
   phase_state_IDs = []
   amp_state_IDs = []
   calibrator_state_IDs = []
   pointing_state_IDs = []
   for state_ID in range(0,len(intents)):
      state_intents = intents[state_ID].rsplit(',')
      for intent in range(0,len(state_intents)):
         scan_intent = state_intents[intent].rsplit('#')[0]
         subscan_intent = state_intents[intent].rsplit('#')[1]
         if (scan_intent == 'CALIBRATE_BANDPASS'):
            bandpass_state_IDs.append(state_ID)
            calibrator_state_IDs.append(state_ID)
         elif (scan_intent == 'CALIBRATE_DELAY'):
            delay_state_IDs.append(state_ID)
            calibrator_state_IDs.append(state_ID)
         elif (scan_intent == 'CALIBRATE_FLUX'):
            flux_state_IDs.append(state_ID)
            calibrator_state_IDs.append(state_ID)
         elif (scan_intent == 'CALIBRATE_POLARIZATION'):
            polarization_state_IDs.append(state_ID)
            calibrator_state_IDs.append(state_ID)
         elif (scan_intent == 'CALIBRATE_AMPLI'):
            amp_state_IDs.append(state_ID)
            calibrator_state_IDs.append(state_ID)
         elif (scan_intent == 'CALIBRATE_PHASE'):
            phase_state_IDs.append(state_ID)
            calibrator_state_IDs.append(state_ID)
         elif (scan_intent == 'CALIBRATE_POINTING'):
            pointing_state_IDs.append(state_ID)
            calibrator_state_IDs.append(state_ID)
