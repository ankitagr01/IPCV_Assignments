Ankit Agrawal (2581532)
Anika Fuchs (2580781)
Akshay Joshi (2581346)

4a) Attached image kodim23s2.ppm, kodim23s4.ppm, kodim23s8.ppm for s=2, s=4 and s=8 respectively.


4b) After reducing the resolution of both the chroma channels, we dont see any significant degrade in the quality of the image. 
And it is expected as, the luma channel Y contains many details and hence we dont subsample it, and the chroma channels contains less details, hence we only subsample them which result in no significant visual deterioration. 
for S=1 its almost same as the original image, visually i cant determine any difference, but as we increase the value of S, we can see very slight change in pixels deteriorating.


4c) Given each channel uses 8 bits, i.e. total of 24 bits are required for an original YCbCr image.
But when we subsample, the Y channel still remains of 8 bits, we only subsample the Cb and Cr channels with the factor of S.
So accordingly, 

For S=2: 
Total bits = 8 + 8/2 + 8/2 = 16 bits

For S=4: 
Total bits = 8 + 8/4 + 8/4 = 12 bits

For S=8: 
Total bits = 8 + 8/8 + 8/8 = 10 bits

