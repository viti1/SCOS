# Optical Neuroimaging Lab Bar-Ilan University
## 1. PlotSCOSvsTime.m
_Calculate Speckle Contrast over Time_  
### Instructions:  
* Add 'func' folder to your matlab path.  
* Run PlotSCOSvsTime.m function. You can run it in a GIU mode or as command with parameters (for more details on command mode see the function header).   
* For GUI : Choose the folder with the recording. If in this folder there is only one .avi file, or a series of .tiff files
it is loaded as a recording. If there are multiple .avi files in the folder, new window is opened and you are 
required to choose the recording file. Otherwise ( no .avi of .tiff files in the folder) error message is presented.
Then, the first frame is imaged and you need to mark a circle on that image where is the region over which calculation
should be performed.    
Next you need to set the window size for culculation - usually 7 or 9 ( better to be an odd number).
That it. Wait few seconds (depends on the recording length) and a figure with SCOS results will appear.
It is automatically saved in the recording folder, as well as a .mat file with the vectors (SCOS vs Time).	

## 2. SCOSloop.m
_run PlotSCOSvsTime over all records in a folder_
I'ts a scrips and available only in command line. 
Change the 'recordsFolder' variable at the beginning of the code to your folder name.
You can change the window size too.

