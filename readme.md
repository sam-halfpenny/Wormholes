This Project Uses Non-Euclidean Raytracing techniques to find null geodesics (the paths light takes) Around the curved spacetime of wormholes described by the Ellis-Bronnikov Wormhole metric.
This code is not efficient in the slightest currently and needs some work to improve it.

To run this project you must change the directories in line 15-23 of wormholes.py to your own. Secondly the code generates a ~1GB lookup table on first run (this will take a long time), every run after that will use the table and so will be much faster.
The benefit of the table is that it will allow rendering from any angle at any radius with any FOV possible in much shorter rendering times (again, further problems need to be addressed in terms of efficiency of this system).

Hope you enjoy! :)