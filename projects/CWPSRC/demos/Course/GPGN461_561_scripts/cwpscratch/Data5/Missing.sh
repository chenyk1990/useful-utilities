#! /bin/sh

# Shell script for dealing with missing traces 

# use: suchart to find missing traces
#

#suchart < seismic.su key1=sx key2=gx |
#xgraph n=120120 linewidth=0 label1="sx" label2="gx" marksize=2 mark=8 &

# by zooming in on the plot, we can find that traces are
# missing between shot positions
# 5187 and 5262
# 5387 and 5487
# 14412 and 14512
# 22162 and 22262


# We will replace the missing shot gathers with the average of the nearest
# neighboring shots gathers

### missing shot gathers to be replaced:
## 5212
## 5237

# capture 5187 and 5262
suwind key=sx min=5187 max=5262 < seismic.su > junk1.su


# sort and stack into an average shot gather
susort dt offset < junk1.su > junk2.su
sustack key=offset < junk2.su > average5187_5262.su

## make shot 5212
#
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=5212,180,0 < average5187_5262.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot5212.su

## make shot 5237
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=5237,181,0 < average5187_5262.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot5237.su


# view concatenation of new shots and seismic.su to see if holes are filled
#cat shot5237.su shot5212.su seismic.su |
#suchart  key1=sx key2=gx  |
#xgraph n=120120 linewidth=0 label1="sx" label2="gx" marksize=2 mark=8 &

# 5387 and 5487
# capture 5387 and 5487
suwind key=sx min=5387 max=5487 < seismic.su > junk1.su

susort dt offset < junk1.su > junk2.su
sustack key=offset < junk2.su > average5387_5487.su

## make shot 5412
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=5412,188,0 < average5387_5487.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot5412.su

## make shot 5437
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=5437,189,0 < average5387_5487.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot5437.su

## make shot 5462
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=5462,190,0 < average5387_5487.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot5462.su

# view concatenation of new shots and seismic.su to see if holes are filled
#cat shot5462.su shot5437.su shot5412.su shot5237.su shot5212.su seismic.su |
#suchart  key1=sx key2=gx  |
#xgraph n=120120 linewidth=0 label1="sx" label2="gx" marksize=2 mark=8 &



# 14412 and 14512
# capture 14412 and 14512
suwind key=sx min=14412 max=14512 < seismic.su > junk1.su
susort dt offset < junk1.su > junk2.su
sustack key=offset < junk2.su > average14412_14512.su

## make shot 14437
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=14437,549,0 < average14412_14512.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot14437.su

## make shot 14462
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=14462,550,0 < average14412_14512.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot14462.su

## make shot 14487
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=14487,551,0 < average14412_14512.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot14487.su

## view concatenation of new shots and seismic.su to see if holes are filled
#cat shot14487.su shot14462.su shot14437.su shot5462.su shot5437.su shot5412.su shot5237.su shot5212.su seismic.su |
#suchart  key1=sx key2=gx  |
#xgraph n=120120 linewidth=0 label1="sx" label2="gx" marksize=2 mark=8 &


# 22162 and 22262
# capture 22162 and 22262
suwind key=sx min=22162 max=22262 < seismic.su > junk1.su


susort dt offset < junk1.su > junk2.su
sustack key=offset < junk2.su > average22162_22262.su

## make shot 22187
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=22187,859,0 < average22162_22262.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot22187.su

## make shot 22212
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=22212,860,0 < average22162_22262.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot22212.su

## make shot 22237
# set sx,ep,nhs header fields fields compute gx from sx and offset
sushw key=sx,ep,nhs a=22237,861,0 < average22162_22262.su | suchw key1=gx key2=sx  key3=offset a=0 b=1 c=1 > shot22237.su


# view concatenation of new shots and seismic.su to see if holes are filled
cat shot22237.su shot22212.su shot22187.su shot14487.su shot14462.su shot14437.su shot5462.su shot5437.su shot5412.su shot5237.su shot5212.su seismic.su | suchart  key1=sx key2=gx  | xgraph n=120120 linewidth=0 label1="sx" label2="gx" marksize=2 mark=8 &

# concatenate, sort to shot gathers, reset cdp field.
cat shot22237.su shot22212.su shot22187.su shot14487.su shot14462.su shot14437.su shot5462.su shot5437.su shot5412.su shot5237.su shot5212.su seismic.su > seismic1.su 

susort < seismic1.su sx offset > seismic2.su

# Rationale of the next paragraph. There are a total of 2142 cdps 
# in the data set. We first set the cdp field to the value of 10 times
# the midpoint, by averaging sx*10 and gx*10. 
# We then shift all of the cdp values by 16060 which makes the value 
# of the first cdp=125, then we divide all cdp values by 125
# which is 10 times the spacing of 12.5 between midpoints
suchw key1=cdp key2=sx key3=gx b=10 c=10 d=2 < seismic2.su |
suchw key1=cdp key2=cdp key3=cdp a=-16060 b=1 c=0 |
suchw key1=cdp key2=cdp key3=cdp a=0 b=1 c=0 d=125  > seismic3.su

# sort into cdps
susort cdp offset < seismic3.su > seis_repaired.cdp.su




exit 0

 


