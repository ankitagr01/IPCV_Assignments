IPCV Assignment 01
T1

Team Members: 
Anika Fuchs (2580781)
Ankit Agrawal (2581532)
Akshay Joshi (2581346)


4b) After integrating uniform distribution noise while quantisation, the quality of the image has been degraded which was expected on addition of noise.
While comparing from 4a for q=3, we can see the image quantised with noise is more degraded than the one quantised without noise. 
It is important to simulate noise as we understand the impact of noise and that can help to figure who to denoise an image. 

4c) Floyd-Steinberg dithering is a error diffusion algorithm which helps to improve the visual quality of the image while quantising. 
Here, we distribute the quantisation error of a pixel to its neighbouring pixels which are yet to be processed.
So when a image is rounded down, it adds the error to its neighbours which are now more likely to be rounded up, hence decreasing the overall quantisation error. 
Therefore we can see that Floyd-Steinberg dithering method actually helps maintaining the quality of the image. 
