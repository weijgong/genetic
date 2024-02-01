import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

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
    barh_param[4].append('y')

    barh_param[0].append(c)
    barh_param[1].append(ows[i][1]-ows[i][0])
    barh_param[2].append(0.8)
    barh_param[3].append(ows[i][0])
    barh_param[4].append('g')

    ticks_y.append(f"")
    ticks_y.append(f"Target no {i+1}")

    c+=2

plt.figure(figsize=(12,4))
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
plt.savefig("/home/gwj/genetic/sat_algorithm/visual_rescheduled_vtw_ows.png")