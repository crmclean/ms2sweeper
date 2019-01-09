# 2019-01-05

This morning in the shower the lightbulb went off. I an use the interpolated parent intensity of the parent ion during the time of the scan was taken or in the time right before the scan was taken as Liz mentioned. 

# 2019-01-04

I got an idea of something that would be cool to add. Since scans that come from the same sample have already been separated, it may be good to use their occourance as a rt threshold. This would assume that the rt deviation is less than the dynamic exclusion window + whatever more time it took for the next scan to appear. 

Thinking about the future, since UPLC is becoming the standard, this might not be such a bad idea. 

# 2018-01-03

The package works and is totally operational. Next I have to do unit testing to get it to upload to bioconductor. 

Otherwise, its probably good that I just wait to get Liz's feedback on the paper I proposed we'd write before moving forward. 

# 2018-01-02

Happy new year R. Today I will finish this package's beta version. I am currently working on the checklist below. I also need to finish up the function to export the spectra to make it ready to run for metfrag. 

## To do

Write an md file for all functions that are going to remain within the code using Roxygen2. 


# 2018-12-21

AWESOME - the hard part is over and it seems to work. At this moment, I have completed the algorithm required to merge ms2 spectra between samples. I still need to do a few things.

COMPLETE: Make probability thresholds input parameters within alignMS2s - specifically within the poisson and uniform distribution parts of the function. 

COMPLETE: Write a function that does extra filtering if needed. Something that only considers peaks above a certain intensity of the max peak or looks for peaks that occur in 2+ scans...

COMPLETE: Third - clean up the sweeper obj by deleting the raw data after going through harvest ms2 function.  

Put together the package and make sure it works on another computer.  

SENT LIZ AN EXAMPLE OF THE TYPE OF PAPER I WILL SEND: Finally - get an outline going to write the paper.  

# 2018-12-18

## State of affiars

1) the group thing didn't work since all features formed a cluster. The new approach probably will require I just check similarity of MS2 spectra collected in a pairwise manner based on increasing rt. 

General overview:
0.5) put the first feature into a group. 
1) check first pair of ms2s
2) if they are the same, join then and check the next one
3) if they are not, create a new group for the second, and check the next one against all groups. 

## updates

Ok - the algorithm is almost done. The functions for grouping scans across samples is all that is left. It has one problem that remains to be solved - confirm that two scans are equal with out doing a MS2 similarity check on each pair of data. 

## current idea

1) separate things into rt groups. 
2) assign communities between things that have similar rts within an allowed error window. Results in an assignment of a family ID.
3) check ms2 spectra for similarities between scans to confirm proper membership. Instead of doing this in an exhaustive manner, check pairs of adjacent features in retention time. (do this for each individual community)
4) reassign features that don't match to distinct communities. 

# 2018-12-10

OK - I need to work on the behavior of the parentPurity function of this algo. There is a disconnect with how I am using the model to predict the behavior. Maybe the play is to find slope between common MS2 features. The key would be to identify IF an MS2 is coming from the same feature. 

## Idea:
1) group elements of the MS2
2) find the slope between the intensity of two adjacent ms2 features
3) check to see which two scans contribute to this by matching their slope

## Thing to consider:
Save this part of the problem for the cleanMS2 function? 
