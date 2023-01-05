# POSSUM_status_monitor
Repository for code relating to monitoring and presenting the status of the POSSUM survey.

#### Jan 5 2022, Cameron:  
I've built an initial prototype of the status monitor, as I envision it.  
In this version, the core is a Google Sheet that holds all the information on the status, for both observations and tiles.
Right now it's only set up for Band 1 (EMU). It has blank sheets for Pilot data, but those are not formatted yet, and
there is no code to automatically update those. (I'll probably just manage those in my testing notebook.)

The Python script in this repository has code for various kinds of updates:
- **Check for new observations.** This reads from the observatory (EMU) status Google Sheet, updates new observations
with SB numbers and observation dates.
- **Update observation with AusSRC processing status.** Flags an observation as having completed AusSRC processing (creating partial tiles in that pipeline).
- **Update tile with AusSRC processing status.** Flags a tile as having completed AusSRC procesing and upload to storage/archive (a fully mosaicked tile).
- **Create plots.** Makes the plots showing the survey coverage/progress, both observations and tiles, in equatorial and Galactic coordinates.
- **Generate an Aladin Lite viewer.** Creates a webpage with an Aladin Lite applet, with overlays showing the observation and tile statuses.

Features that are still missing:
- *Everything for band 2 (WALLABY).* Waiting on access to WALLABY observatory status sheet, as well as full field list (with all beam positions).
- *Validation status updates.* Need to figure out if there's a way to auto-query validation status. Checking if cubes are public?
- *Pipeline status updates.* Need to write pipeline module(s) that update the tile status through the pipelines.
- *Confirm things work with re-observations.* If a field is reobserved, will that be captured in the spreadsheet update? Need to check.
- *Pilot data.* Will probably want sheets for band 1 and 2 separately.
- *Website* to hold figures and Aladin lite page. HTML page to hold figures?


#### Links:  
Spreadsheet for automatic update: https://docs.google.com/spreadsheets/d/1sWCtxSSzTwjYjhxr1_KVLWG2AnrHwSJf_RWQow7wbH0/edit#gid=0

Manual spreadsheet: https://docs.google.com/spreadsheets/d/1D2ZoaWD9jFVZBf8njYLMhoy0vyi5EuT2e37SwYvHQ8E/edit#gid=1163143711


EMU Observatory status sheet: https://docs.google.com/spreadsheets/d/1HjSWDvwknndCi5PxTP4aDHYAdQiSliRYFVwZGf8YZJw/edit#gid=0

EMU Status page: 
- http://www.emu-survey.org/progress/
- http://www.emu-survey.org/progress/aladin.html
- http://www.emu-survey.org/progress/table.html


