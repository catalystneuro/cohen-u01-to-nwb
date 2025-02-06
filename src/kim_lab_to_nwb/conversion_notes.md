

## Structure of the matlab file


# data(1,:) is time
# data(2,:) is left wingbeat
# data(3,:) is left-right wingbeat
# data(4,:) is x-position of the visual pattern
# data(5,:) is y-position of the visual pattern
# data(6,:) is 2-photon frame synchronization signal (1 pulse corresponds to 1 frame)
# data(7,:) is behavior camera signal (1 pulse corresponds to 1 frame)
# data(8,:) indicates the start of a stimulus (it is empty in this example)

# Questions for the meeting on 2024-11-06
* How is the data in the matlab aquired? which DAQ?
The DAQ is the following:
[Link to device](https://www.digikey.com/en/products/detail/ni/782258-01/12817857)
* The data that you shared with us? is this associated with a paper? is this something you are interested on sharing on dandi? What are the expectations of the code 
No paper, what they want is some flexible set of routines that they can upload to DANDI.
* The tiff format, what kind of microscope is it? Scanimage? If you could send me the model that would be useful.
The microscope is self-made but for data acquisition they use Scanimage. They will send some metadata.
* How is the stimuli presented? That is, I have two images, are they presented on time?`
Those are images, a new dataset with the stimuli will be send to us.
* How to synchronize the Videos, we should be aligned, we are dtecting crossin a threshold for hits.
They are synchronized using the data in the DAQ.
* Are the timestamps of the fluoresence traces and the stimuli the ones that are in the DAQ as well?
I see that all of them have the same number of samples.
* What are the emission oand excitation lambda?

The data from the stimuli (in the raw folder `Pattern_1_stripe_24.mat`) looks like this:

```
x_num 161
y_num 14
num_panels 24
row_compression 0
gs_val 1
Pats shape=(24, 96, 161, 14)
Panel_map shape=(2, 12)
data shape=(432768,)
```

It also has a BitMaxIndex that looks more convoluted:

```
BitMapIndex {'Panel_ID': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], 'row_range': [array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8)], 'column_range': [array([73, 74, 75, 76, 77, 78, 79, 80], dtype=uint8), array([49, 50, 51, 52, 53, 54, 55, 56], dtype=uint8), array([25, 26, 27, 28, 29, 30, 31, 32], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([65, 66, 67, 68, 69, 70, 71, 72], dtype=uint8), array([41, 42, 43, 44, 45, 46, 47, 48], dtype=uint8), array([17, 18, 19, 20, 21, 22, 23, 24], dtype=uint8), array([89, 90, 91, 92, 93, 94, 95, 96], dtype=uint8), array([57, 58, 59, 60, 61, 62, 63, 64], dtype=uint8), array([33, 34, 35, 36, 37, 38, 39, 40], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([81, 82, 83, 84, 85, 86, 87, 88], dtype=uint8), array([73, 74, 75, 76, 77, 78, 79, 80], dtype=uint8), array([49, 50, 51, 52, 53, 54, 55, 56], dtype=uint8), array([25, 26, 27, 28, 29, 30, 31, 32], dtype=uint8), array([1, 2, 3, 4, 5, 6, 7, 8], dtype=uint8), array([65, 66, 67, 68, 69, 70, 71, 72], dtype=uint8), array([41, 42, 43, 44, 45, 46, 47, 48], dtype=uint8), array([17, 18, 19, 20, 21, 22, 23, 24], dtype=uint8), array([89, 90, 91, 92, 93, 94, 95, 96], dtype=uint8), array([57, 58, 59, 60, 61, 62, 63, 64], dtype=uint8), array([33, 34, 35, 36, 37, 38, 39, 40], dtype=uint8), array([ 9, 10, 11, 12, 13, 14, 15, 16], dtype=uint8), array([81, 82, 83, 84, 85, 86, 87, 88], dtype=uint8)]}
```

## For the ROI there are three fields

This is in the analysis folder (`ROI_20240108b_00003.mat`), the fields are:

```

dict_keys(['x_cor', 'reference_image', 'y_cor', 'total_ROI'])
```

in the data that was shared with us their shapes are:

```
total_ROI: 17
x_cor: shape =(17, 3)
y_cor: shape=(17, 3)
```

What are the three things, vertices? So the ROIs are triangles.


## Imaging data
This is scan image data with a self-made microscope. The specification of the tiff file metadata can be found here:

https://docs.scanimage.org/Appendix/ScanImage%2BBigTiff%2BSpecification.html#scanimage-bigtiff-specification

