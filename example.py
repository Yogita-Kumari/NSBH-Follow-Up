#To run it 
# python example.py   --output-path /home/yogitakumari/Documents/NSBH_lc/gwemopt_runs/example_output   --numDays 1   --coinc-time 1264316116   --skymap /home/yogitakumari/Documents/NSBH_lc/gwemopt_runs/100.fits --lightcurveFiles /home/yogitakumari/Documents/NSBH_lc/gwemopt_runs/lightcurve_mejdyn=0.01mej_wind=0.01theta=45.dat

#Importing libraries
import os,sys
import tempfile
from pathlib import Path
import glob
import numpy as np
import pandas as pd
from astropy.time import Time
from gwemopt.io.schedule import read_schedule
from gwemopt.run import run
import logging 
import subprocess
import argparse


np.random.seed(1235)


parent_dir = Path(__file__).parent.absolute()

logging.basicConfig(format='%(asctime)s : %(message)s', level=logging.DEBUG)
parser = argparse.ArgumentParser()
# parser.add_argument("--verbose", action="count")
parser.add_argument("--output-path",
                    help="Path to the output file")
parser.add_argument("--numDays",
                    help="Number of days of observation")
# parser.add_argument("--Ntiles",
#                     help="Number of tiles used") #Yogita: No need in case of moc
parser.add_argument("--coinc-time",
                    help="Input the coincidence time from the fits")
parser.add_argument("--skymap",
                    help="fits skymap file ")

parser.add_argument("-l", "--lightcurveFiles",
                    help="lightcurve .dat file")


args = parser.parse_args()

output_dir = args.output_path
coinc_time = args.coinc_time
skymap = args.skymap

# gwemopt settings 
numDays = args.numDays
# ntiles_used = args.Ntiles

lightcurve=args.lightcurveFiles

logging.info('Obtaining fits %s', skymap)

results_dir = Path(output_dir)



# telescope settings  
telescope_list = [
    (
        "ZTF",
        [
            # "--doCoverage",
            # "--coverageFiles",
            # os.path.join(results_dir, "coverage_ZTF.dat"), #These are of use only if --doSchedule is not there
            "--filter", "r",
            "--exposuretimes","30",
            "--doChipGaps",
        ],
     ),
        ( "WINTER",
         [
            #  "--doCoverage",
            #  "--coverageFiles",
            #  os.path.join(results_dir, "coverage_WINTER.dat"),
             "--filter", "J",
             "--exposuretimes","450",
         ],

        ),
        ( "LSST",
            [
                # "--doCoverage",
                # "--coverageFiles",
                # os.path.join(results_dir, "coverage_LSST.dat"),
                "--filter", "r",
                "--exposuretimes","30",
            ],

        )
]

#         '--Ntiles',str(ntiles_used),

for telescope, extra in telescope_list:
    output_dir_telescope = results_dir.joinpath(telescope)
    args = [
        f"-t",
        telescope,
        "-o",
        str(output_dir_telescope),
        "-e",
        str(skymap), 
        # "--doSkymap", #byme
        "--doTiles",
        "--plots", "tiles", #, "skymap","coverage", "efficiency" , #Added
        #"--doPlots", #byme
        "--doSchedule",
        "--timeallocationType","powerlaw",
        #"--powerlaw_cl", "0.9", #byme
        "--doSingleExposure",
        "--doAlternatingFilters",
        "--doBalanceExposure",
        "--geometry", "3d",
        # '--do3D',   #byme
        '--gpstime', str(Time(coinc_time, format='gps', scale='tcg')),
        "--tilesType","moc",
        '--Tobs', '0.0,' + str(numDays),
        "--doReferences",
        # "--doObservability",
        "--doOrderByObservability", #Added
        "--doEfficiency",
        "--modelType", "file", #Added
        "--Ninj","1",  #Added
    ] + extra

    logging.info('Running gwmeopt ....')

    try:
        #byme next 3 lines
        # Run the command and capture both stdout and stderr
        # process = subprocess.Popen(run(args), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        # stdout, stderr = process.communicate()
        logging.info(f"Command being run: {args}")

        run(args)
        
        #byme next 4 lines
        # # Log the command output (stdout and stderr)
        # logging.info('Command output (stdout):\n%s', stdout)
        # if stderr:
        #     logging.warning('Command output (stderr):\n%s', stderr)


        check_files = [
        f"schedule_{telescope}.dat",
#        f"tiles_coverage_int_{telescope}.txt", #byme
#         f"coverage_{telescope}.dat",
    ]
       
        logging.info('gwemopt run completed..')
        logging.info('gwemopt testing begins..')
        
        for i, file_name in enumerate(check_files):
            logging.info(f"Testing: {file_name}")
            
            
            new_schedule = read_schedule(Path(output_dir_telescope).joinpath(file_name))
            expected_schedule = read_schedule(
                output_dir_telescope.joinpath(file_name)
            )

            pd.testing.assert_frame_equal(
                new_schedule.reset_index(drop=True),
                expected_schedule.reset_index(drop=True),
            )

        
        # Test the extra efficiency/coverage files
        extra_test_files = [
#            "map.dat",  #byme
            "summary.dat",
#             "coverage_{telescope}.dat",
        ]

        for extra_test_file in extra_test_files:
            logging.info(f"Testing: {extra_test_file}")
            new = pd.read_table(
                Path(output_dir_telescope).joinpath(extra_test_file), sep=r"\s+"
            )
            expected = pd.read_table(
                output_dir_telescope.joinpath(extra_test_file), sep=r"\s+"
            )
            pd.testing.assert_frame_equal(
                new.reset_index(drop=True),
                expected.reset_index(drop=True),
            )

        logging.info('gwemopt run completed..')


    except Exception as e:
        logging.error('Error running the command: %s', str(e))