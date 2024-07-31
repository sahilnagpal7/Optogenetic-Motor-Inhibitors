# Optogenetic-Motor-Inhibitors
Codes used for data analysis in Nagpal et al., 2024

# Instructions for the velocity code:

CSV file generated from ImageJ plugin - TrackMate (code compatible with CSV generated from the older version of TrackMate) is used in the first code. Manually add x and y coordinates of the cell center to V1 and W1 in the CSV file. 

First code - "Velocity_analysis_part1.m":
Things to do before running the script:
Mention the directory of your CSV data - line 23
Mention the name of the CSV file that you want to analyze - line 25
Things to do when running the script:
Once you run this code, do the following:
For this question, "Do you want to manually time series for the conditions? 1 for yes, 2 for no : ", answer 1
Enter the number of conditions you want!So you will have to mention 4 time series values. Type 0 (lower), 120 (upper), 120 (lower) and then 255 (upper).
The code will process the velocity and sliding window time. You will have to save the "velocity_1" and "Slid_time" manually from the Workspace! After using this code, proceed to "Velocity_analysis_part2.m" code to plot the velocity time series plot.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Second code - "Velocity_analysis_part2.m":
Things to do before running the script:
Mention the directory of your CSV data - line 10
You can play with different values of average velocity time window (30, 20, 10 or 5 - the time window before the inhibition begins). You can change the window by mentioning different values in this part of the code line (bolded number - 85 in this case): "condition_indx = find(Sliding_time_processed{tn} > 85 & Sliding_time_processed{tn} <= 115);. This is in line 50.
You can have different values of the neutral velocity range like +/- 0.01, 0.02 or 0.03 (0.04 is not recommended cause the max positive velocity is around 0.045 I think). You can change this range in the line 71 (positive velocity) and line 84 (negative velocity).
Things to do when running the script:
Once you run this code, do the following:
The code asks you the name of the velocity file you saved manually. Write the name of the file (please don't write the extension)
Next, the code asks you, "How many conditions do you want to consider?". You can have 1 condition (Dark1 alone; 2 values to mention - 0 and 120), 2 conditions (Dark1 and Lit; 4 values to mention - 0, 120, 120 and 255
