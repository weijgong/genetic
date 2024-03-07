from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# create new figure, axes instances.
pos_list = []
with open("/home/gwj/genetic/sat_data/Target1_LLA_Position.csv") as f:
    pos_list = f.readlines()
pos_list = [i.strip().split(',')[:-1] for i in pos_list]

fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])

bound_los,bound_loe,bound_las,bound_lae = 20.,35.,0,15.
m = Basemap(llcrnrlon=bound_los,llcrnrlat=bound_las,urcrnrlon=bound_loe,urcrnrlat=bound_lae,
            rsphere=(6378137.00,6356752.3142),
            resolution='l',projection='merc',
            lat_0=(bound_las+bound_lae)/2,lon_0=(bound_los+bound_loe)/2)

maped_list = []

for x,y in pos_list:
    x,y = float(x),float(y)
    # xpt: longitude ypt: latitude
    xpt,ypt = m(y,x)
    maped_list.append([xpt,ypt])
    # 绘制各个目标在地图上的点位置
    m.plot(xpt,ypt,'k.')

for i in range(len(maped_list)-1):
    # 计算箭头的中点，用于定位箭头标签  
    mid_x = (maped_list[i][0] + maped_list[i+1][0]) / 2  
    mid_y = (maped_list[i][1] + maped_list[i+1][1]) / 2  
    
    # 使用annotate绘制箭头,表示观测的具体顺序
    plt.annotate('', xy=(maped_list[i+1][0], maped_list[i+1][1]), 
                xytext=(maped_list[i][0], maped_list[i][1]),  
                xycoords='data', textcoords='data',  
                arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3,rad=.6',  
                                lw=1, color='black')
                ) 

m.drawcoastlines()
m.fillcontinents()
# labels = [left,right,top,bottom]
# 绘制竖直方向的线（经线）
m.drawmeridians(np.arange(bound_los, bound_loe, 2.5),labels=[False,False,False,True])
# 绘制水平方向的线（纬线）
m.drawparallels(np.arange(bound_las, bound_lae, 2.5),labels=[True,True,False,False])

ax.set_title('The sequence of observe target')
plt.savefig("visual_observe_sequce_src.png")