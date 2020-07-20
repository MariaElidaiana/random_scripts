from glob import glob
import pandas as pd

#Get the list of logfiles
TMlogdir  = '/data/des41.a/data/desgw_public/TreasureMap/logs/'
logdir    = glob(TMlogdir + "*/")
filenames = [glob(i + "*.log")[0] for i in logdir]

#Get the info from the logfiles
user, event, date, bands = [],[],[],[]
for i in range(len(filenames)):
    f = open(filenames[i])
    lines = f.readlines()
    user_event_str = lines[3].split(" ")
    bands_str = lines[11].split(" ")
    user.append(user_event_str[5].replace("[","").replace("]",""))
    event.append(user_event_str[9].replace("\n", "") )
    date.append(user_event_str[1])
    bands.append(bands_str[9] )

#Save as csv
df = pd.DataFrame({"user": user, "event": event, "date": date, "bands":bands}, columns=['user', 'event', 'date', 'bands'])
df.to_csv("./DESGW_TreasureMap_submitted.csv", sep=',',index=False)
