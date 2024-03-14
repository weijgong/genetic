import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

lines = []
vtws = []
ows = []
mode = []
with open("/home/gwj/genetic/sat_algorithm/time_windows_data.dat",'r') as f:
    lines = f.readlines()

oceanic_colors = ['#8ED6FF', '#5D9CEC', '#265A84', '#48CFAD', '#16A085']
# print(lines)
for i in lines:
    tmp = i.strip('\n').split(',')
    tmp = [float(j) for j in tmp]
    vtws.append([tmp[0],tmp[1]])
    ows.append([tmp[2],tmp[3]])
    mode.append(tmp[4])
    # print(tmp)
# print(vtws,ows,mode)
# barh_param: y,width,height,left
barh_param = [[],[],[],[],[]]
ticks_y = []
c = 0
for i in range(len(vtws)):        
    barh_param[0].append(c-0.8)
    barh_param[1].append(vtws[i][1]-vtws[i][0])
    barh_param[2].append(0.8)
    barh_param[3].append(vtws[i][0])
    barh_param[4].append(oceanic_colors[1])

    barh_param[0].append(c)
    barh_param[1].append(ows[i][1]-ows[i][0])
    barh_param[2].append(0.8)
    barh_param[3].append(ows[i][0])
    barh_param[4].append(oceanic_colors[0])

    ticks_y.append(f"")
    ticks_y.append(f"Target no {i+1}")

    c+=2

plt.figure(figsize=(26,10))
plt.subplot(2,1,1)
plt.barh(
    y=barh_param[0],
    width=barh_param[1],
    height=barh_param[2],
    left=barh_param[3],
    align='center',
    color = barh_param[4],
    label = ['visual windows','observe windows']+["" for i in range(len(barh_param[0])-2)]
    )

lines = []
free_windows = []
with open("/home/gwj/genetic/sat_algorithm/free_time_windows_data.dat",'r') as f:
    lines = f.readlines()
for i in lines:
    tmp = i.strip('\n').split(',')
    tmp = [float(j) for j in tmp]
    free_windows.append(tmp)

cur_left = 0
for i in range(len(free_windows)):
    print(free_windows[i])
    left_  = free_windows[i][0]
    width_ = free_windows[i][1] - free_windows[i][0]
    if i == 0:
        label_ = ['free observation windows']
    else:
        label_ = ['']
    plt.barh(
    y=barh_param[0][-1]+1.5,
    width=width_,
    height=0.8,
    left=left_,
    align='center',
    color =oceanic_colors[2],
    label = label_
    )

plt.grid()
plt.legend(prop={'size': 6})
plt.title("satellite schedule")
plt.xlabel("Timeline")
plt.ylabel("target no")
plt.yticks(barh_param[0]+[barh_param[0][-1]+1.5],labels=ticks_y+['Free windows'])
# plt.savefig("/home/gwj/genetic/sat_algorithm/visual_vtw_ows.png")



lines = []
vtws = []
ows = []
mode = []
with open("/home/gwj/genetic/sat_algorithm/rescheduled_time_windows_data.dat",'r') as f:
    lines = f.readlines()

# print(lines)
for i in lines:
    tmp = i.strip('\n').split(',')
    tmp = [float(j) for j in tmp]
    vtws.append([tmp[0],tmp[1]])
    ows.append([tmp[2],tmp[3]])
    mode.append(tmp[4])
    # print(tmp)
# print(vtws,ows,mode)
# barh_param: y,width,height,left
barh_param = [[],[],[],[],[]]
ticks_y = []
c = 0
for i in range(len(vtws)):        
    barh_param[0].append(c-0.8)
    barh_param[1].append(vtws[i][1]-vtws[i][0])
    barh_param[2].append(0.8)
    barh_param[3].append(vtws[i][0])
    barh_param[4].append(oceanic_colors[1])

    barh_param[0].append(c)
    barh_param[1].append(ows[i][1]-ows[i][0])
    barh_param[2].append(0.8)
    barh_param[3].append(ows[i][0])
    barh_param[4].append(oceanic_colors[0])

    ticks_y.append(f"")
    ticks_y.append(f"Target no {i+1}")

    c+=2

# plt.figure(figsize=(12,4))

plt.subplot(2,1,2)
plt.barh(
    y=barh_param[0],
    width=barh_param[1],
    height=barh_param[2],
    left=barh_param[3],
    align='center',
    color = barh_param[4],
    label = ['visual windows','observe windows']+["" for i in range(len(barh_param[0])-2)]
    )

plt.grid()
plt.legend(prop={'size': 6})
plt.title("rescheduled satellite timeline digram")
plt.xlabel("Timeline")
plt.ylabel("target no")
plt.yticks(barh_param[0]+[barh_param[0][-1]+1.5],labels=ticks_y+['Free windows'])



plt.savefig("/home/gwj/genetic/sat_algorithm/visual_all_vtw_ows.png")