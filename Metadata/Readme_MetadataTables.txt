'LoggerTables.xlsx' -- gives information about each raw data table from logger

1. "Site": Location name for data collection
2. "Date": Date for data collection (yyyy-mm-dd)
3. "Logger": Logger name for data collection
4. "Table": Data table name for data collection
5. "Variables": Comma-separated list of variable names for data table.  Guide to variable names is in "LoggerVariables.xlsx"
6. "ErrorCode": Set to 1 if timing errors may be present in table.  Otherwise set to 0.

---

'LoggerTimes.xlsx' -- gives information about time intervals for each logger for each day

1. "Site": Location name for data collection
2. "Date": Date for data collection (yyyy-mm-dd)
3. "Logger": Logger name for data collection
4. "StartTime": Beginning of time interval for usable data from logger (HH:MM)
5. "EndTime": End of time interval for usable data from logger (HH:MM)
6. "Delay_Seconds": Duration by which to shift times to account for logger sync.  Positive values are shifted to later times, negative values are shifted to earlier times.  Generally determined based on slave pulses received by loggers.

---

'InstrumentVariables.xlsx' -- gives information instruments and associated variables

1. "InstrumentType": Gives type of instrument, which determines method of processing data
2. "Instrument": Name of instrument
3. "VarNameSpecific": Specific variable name, which ties variable to a specific instrument
4. "VarNameGeneral": General variable name, which is not instrument specific
5. "DataType": Data type for variable values.  "float" for decimal values and "int" for integer values.  This determines method of parsing when importing data from logger files
6. "Units": Units for this variable.  Included in metadata when parsing variable.
7. "Calibration": Equals 1 when no calibration is required.  If value differs from one, then multiply by this value when performing calibration.

---

'InstrumentMetadata.xlsx' -- Gives information about deployment intervals for each day for each instrument at site.  Intervals are only those times with valid data (i.e., no interference with instruments)

1. "Site": Location name for data collection
2. "Date": Date for data collection (yyyy-mm-dd)
3. "StartTime": Beginning of time interval for valid data from instrument (HH:MM)
4. "EndTime": End of time interval for valid data from instrument (HH:MM)
5. "InstrumentType": Gives type of instrument, which determines method of processing data
6. "Instrument": Gives name of instrument
7. "Variables": Gives names of variables for this instrument
8. "StartHeight_m": Height of instrument measured at beginning of day, in meters
9. "EndHeight_m": Height of instrument measured at end of day, in meters
10. "HeightErr_m": Estimated error (+/-) in instrument height measurements, in meters
11. "HeightRef": Set to 0 if height is given relative to ground.  If set to 'L1' or 'L2', then height is given relative to named distance sensor
12. "Longitudinal_m": Downwind displacement of instrument (in m), relative to setup orientation
13. "Spanwise_m": spanwise displacement of instrument (in m), relative to setup orientation
14. "AngleErr_deg": angle of mismatch (in deg) between orientation of setup and actual wind.  Positive angle corresponds to angular shift of wind toward positive y according to right-hand rule.  Can use this value to adjust longitudinal and spanwise positions.
15. "ErrorCode":  Set to 1 if a major problem is detected with instrument during time interval, otherwise set to 0.  This argument allows easy toggling to include / exclude instrument as analysis progresses.

---

'WeightBSNE.xlsx' -- Gives information about BSNE collections

1. "Site" = Location name for data collection
2. "Date": Date for data collection (yyyy-mm-dd)
3. "StartTime": Time when trap was opened (HH:MM)
4. "EndTime": Time when trap was closed (HH:MM)
5. "NameBSNE": Name assigned to BSNE within array
6. "Weight_g": Weight of sample (g)
7. "WeightErr_g": Estimated error associated with weight (g)
8. "BottomHeight_cm": Height from ground to bottom of opening (cm)
9. "BottomHeightErr_cm": = Uncertainty in BSNE heights related to measurement over rough bed and changes in height with time = (+/- 0.5 cm)
10. "HeightBSNE_cm": Height of BSNE opening (cm).  5 cm for standard BSNE, 1 cm for modified BSNE.
11. "WidthBSNE_cm": Width of BSNE opening (cm). 2 cm for both kinds of BSNEs.
12. "Longitudinal_m": Downwind distance of BSNE relative to Wenglor / sonic colocation point (m)
13. "Spanwise_m": Spanwise distance of BSNE direction relative to Wenglor / sonic colocation point (m)
14. "ErrorCode": = 0 if no problems with data.  = 1 if problems suspected.  Assign ErrorCode = 1 to all samples with mass less than 0.1 g.
15. "StartTime_GrainSize": Time for start of grain size analysis for trap (HH:MM). Set to 00:00 if no grain size is available (sample too small, probably).  Time interval may combine multiple BSNE trap times so that sufficient sand available for grain size analysis.
16. "EndTime_GrainSize: Time for end of grain size analysis for trap (HH:MM).  Set to 00:00 if no grain size is available (sample too small, probably).  Time interval may combine multiple BSNE trap times so that sufficient sand available for grain size analysis.