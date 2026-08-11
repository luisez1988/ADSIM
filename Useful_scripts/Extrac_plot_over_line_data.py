#### EXPORT PLOT OVER LINE FOR ALL TIME STEPS ####

from paraview.simple import *
import os

# ------------ USER SETTINGS ------------
# needs existing active plot over line object
plotLine = GetActiveSource()

output_folder = r"C:\Users\zamcr\Documents\GitHub\Anura3D_OpenSource\ADSIM\src\output"    # change to your folder
filename_pattern = "t{step:02d}.csv"   # custom naming pattern

# Create folder if needed
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ------------ GET ACTIVE DATA ------------
source = plotLine
if source is None:
    raise RuntimeError("No active source. Select your data first.")

# ------------ TIME MANAGEMENT ------------
animationScene = GetAnimationScene()
timeKeeper = GetTimeKeeper()

# Force animation scene to match available times
animationScene.UpdateAnimationUsingDataTimeSteps()

# List of all time values
time_steps = timeKeeper.TimestepValues

print("Found", len(time_steps), "time steps.")

# ------------ EXPORT LOOP ------------
for i, t in enumerate(time_steps):
    print(f"Processing time step {i} (t={t})...")

    # Set animation time
    animationScene.TimeKeeper.Time = t

    # Update filter
    plotLine.UpdatePipeline(t)

    # Build file name
    file_name = filename_pattern.format(step=i)
    file_path = os.path.join(output_folder, file_name)

    # Export CSV
    SaveData(file_path, proxy=plotLine)
    print(" -> saved:", file_path)

print("Done!")
