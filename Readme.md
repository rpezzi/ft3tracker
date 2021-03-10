# ALICE 3 Forward Tracker

### 1. Simulation

`$ o2-sim -m FT3 -e TGeant3 -g pythia8 -n 1000`

Output: `o2sim_HitsFT3.root`

### 2. Tracking
`$ root.exe -q ft3Tracker.C+`

Output: `ft3tracks.root`

### 3. Assessment histograms
`$ root.exe -q -b FT3TrackerChecker.C+`

Output: `Fittercheck_ft3tracks.root`
