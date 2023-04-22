#Case-oriented spatial and temporal density scanning
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta
import time


#函数，给定一个病例，统计以该病例为中心的时空立方体内的所有病例。   返回值是该病例为中心的的时空立方体内的所有病例下标数组，第[0]个存的是个数。
def getSTcube(CaseData,case_index,Time_interval:float,Lon_interval:float,Lat_interval:float,case_in_cube_numMax = 10000,Date_FieldName="发病日期",Lon_FieldName="lon",Lat_FieldName="lat"):
    #CaseData 病例表，DataFrame格式，包含所有需要的信息。 必须是按照时间顺序排列的
    #case_index 中心病例在数组中的下标
    #Time_interval 单位为天
    #Lon_interval, Lon_interval 单位为度
    Date_of_illness=CaseData[Date_FieldName];
    Lon_of_illness=CaseData[Lon_FieldName];
    Lat_of_illness=CaseData[Lat_FieldName];
    
    case_in_cube = np.zeros(case_in_cube_numMax+1,dtype=np.int32);
    #case_in_cube[0]为保留空间，用来记录该方格中病例的个数
    #下面的算法会保证，case_in_cube[1]记录的是给定中心病例的坐标

    Case_date=datetime.strptime(Date_of_illness[case_index],"%Y/%m/%d");
    #以所给病例为中心，向两侧统计，直到超出时间阈值范围
    #1.向前搜索：
    for i in range(case_index,-1,-1):
        #for i in range(case_index,-1,-1)，第一个数是'case_index', 保证了 case_in_cube[1]记录的是给定中心病例的坐标
        date=datetime.strptime(Date_of_illness[i],"%Y/%m/%d");
         #print((date-Case_date).days)   #实际上该属性为整数，存在舍入
        if abs(Case_date-date).days <= Time_interval/2:
            #时间上，位于所给病例为中心的时空立方体内
            if abs(Lon_of_illness[i] - Lon_of_illness[case_index]) <= Lon_interval/2:
                if abs(Lat_of_illness[i] - Lat_of_illness[case_index]) <= Lat_interval/2:
                    #意味着此病例位于所给病例为中心的时空立方体内
                    case_in_cube[0] = case_in_cube[0] +1;
                    if case_in_cube[0] > case_in_cube_numMax:
                        print(f'\n\n\nError: the count of cases in one cube(index:{case_index}) is more than case_in_cube_numMax = {case_in_cube_numMax}, please give a bigger case_in_cube_numMax\n\n\n')
                        return case_in_cube;
                    case_in_cube[case_in_cube[0]] = i;

        else:
            #时间上，已经超出所给病例为中心的时空立方体，因为有序，其余病例也将都超过，所以直接退出循环
            break;


    #2.向后搜索：==========================
    for i in range(case_index+1,CaseData.shape[0]):
        date=datetime.strptime(Date_of_illness[i],"%Y/%m/%d");
        
        if abs(date-Case_date).days <= Time_interval/2:
            #时间上，位于所给病例为中心的时空立方体内
            if abs(Lon_of_illness[i] - Lon_of_illness[case_index]) <= Lon_interval/2:
                if abs(Lat_of_illness[i] - Lat_of_illness[case_index]) <= Lat_interval/2:
                    #意味着此病例位于所给病例为中心的时空立方体内
                    case_in_cube[0] = case_in_cube[0] +1;
                    if case_in_cube[0] > case_in_cube_numMax:
                        print(f'\n\n\nError: the count of cases in one cube(index:{case_index}) is more than case_in_cube_numMax = {case_in_cube_numMax}, please give a bigger case_in_cube_numMax\n\n\n')
                        return case_in_cube;
                    case_in_cube[case_in_cube[0]] = i;

        else:
            #时间上，已经超出所给病例为中心的时空立方体，因为有序，其余病例也将都超过，所以直接退出循环
            break;
                
    return case_in_cube;

  
  
  
def Get_ST_Scan(CaseData,Time_interval, Lon_interval, Lat_interval, outbreaks_num_Threshold, 
                Time_Threshold, Lon_Threshold, Lat_Threshold, Date_FieldName, Lon_FieldName, Lat_FieldName, outbreaks_num_Threshold_big = 20,
                case_in_cube_numMax = 1000, case_in_cluster_numMax = 10000,outbteaks_csae_cube_numMax = 10000):
    
    print(f'\n\n\nTime_interval: {Time_interval};  Lon_interval: {Lon_interval};  Lat_interval: {Lat_interval};  outbreaks_num_Threshold: {outbreaks_num_Threshold};  ');
    print(f'Time_Threshold: {Time_Threshold}倍Time_interval;  Lon_Threshold: {Lon_Threshold}倍Lon_interval;  Lat_Threshold: {Lat_Threshold}倍Lat_interval; ');
    print(f'Date_FieldName: {Date_FieldName};  Lon_FieldName: {Lon_FieldName};  Lat_FieldName: {Lat_FieldName};   \n');



    print("\n开始 病例时空密度统计......");
    Process_start_time=time.perf_counter();
    Process_Part_start_time=time.perf_counter();


    #Part1 逐病例为中心进行时空立方体统计
    outbreaks_case_num = 0;
    outbteaks_csae_cube_Case_Indexs=np.zeros((outbteaks_csae_cube_numMax,case_in_cube_numMax+1),dtype=np.int32);  #(num,Index1,Index2,...)  
    #用来记录每个病例为中心的方格里所有病例的Index，通过 ID_of_illness[Index]等即可获取到病例ID等信息。
    #outbteaks_csae_cube_Case_Indexs[i][0]为保留空间，用来记录该方格中病例的个数; 统计函数会保证outbteaks_csae_cube_Case_Indexs[i][1]记录的是中心病例的坐标


    for i in range(CaseData.shape[0]):
        if 0 == i%1000:
            print(f'{time.strftime("%Y/%m/%d %H:%M:%S",time.localtime())}, "Date_of_Case:"{CaseData[Date_FieldName][i]}, 总计{CaseData.shape[0]}条病例,    已经处理 {i} 条病例，标识出{outbreaks_case_num}个超出阈值的以病例为中心的时空立方体......');
        case_cube = getSTcube(CaseData,i,Time_interval,Lon_interval,Lat_interval,case_in_cube_numMax = case_in_cube_numMax,Date_FieldName = Date_FieldName, Lon_FieldName = Lon_FieldName, Lat_FieldName = Lat_FieldName);
        if case_cube[0] > case_in_cube_numMax:
            print(f'\n\n\nError: 第{i}个病例为中心的时空立方体中，病例数为{case_cube[0]}, 超过设定的阈值case_in_cube_numMax: {case_in_cube_numMax},统计无法继续进行。\n请为函数传入一个更大的 case_in_cube_numMax\n\n\n');
            return;
        if case_cube[0] >= outbreaks_num_Threshold:
            outbteaks_csae_cube_Case_Indexs[outbreaks_case_num][0:case_cube[0]+1] = case_cube[0:case_cube[0]+1];
            outbreaks_case_num += 1;
            if outbreaks_case_num >= outbteaks_csae_cube_numMax:
                print(f'\n\n\nError: 超出阈值的以病例为中心的时空立方体个数超过设定的阈值outbteaks_csae_cube_numMax: {outbteaks_csae_cube_numMax},统计无法继续进行。\n请为函数传入一个更大的 outbteaks_csae_cube_numMax\n\n\n');
                return;
    print(f'outbreaks_case_num: {outbreaks_case_num}');

    print(f"完成 面向病例记录的以病例为中心的时空立方体 统计，该部分耗时 {(time.perf_counter() - Process_Part_start_time)/60} 分钟......");
    print(f"总耗时 {(time.perf_counter() - Process_start_time)/60} 分钟......\n\n\n");
    Process_Part_start_time=time.perf_counter();



    
    
    #合并超出阈值的 以病例为中心的时空立方体。   合并算法借鉴了数字图形处理中，标记连接成分的算法。
    print("\n\n\n\n\n开始 合并超出阈值的 以病例为中心的时空立方体......");


    #合并的根本原因：归并一次爆发事件，因爆发时间和空间上的范围不定，有可能一次相关联的聚集爆发跨越几个时空立方体
    #由于没有病原体的基因序列数据，无法从基因上判断不同病例是否为同一病原体导致的传播，只能采取时空临近的原则。
    #合并的准则：基于时空临近原则，合并的准则本质上就是两个立方体小于一定距离。  
    #           采取相邻合并的话，就是意味着时间上一个单位长度、空间上一个单位长度。
    #           其实思维可以更开放一点，把合并的时间距离和空间距离都放宽一点，例如两个单位长度。
    #           并且，没有必要时间上和空间上采取同样倍数的空间长度。  完全可以时间上2倍单位长度，空间上1倍单位长度。


    #合并 超出阈值的 以病例为中心的时空立方体大步骤：
    #1.先合并中心病例间 时空距离 小于合并阈值要求的 以病例为中心的时空立方体
    #2.完成后，再合并有相同病例的 超出阈值的 以病例为中心的聚集爆发事件（当时空距离阈值大于等于1单位时，理论上第一种合并后，不会出现这种合并）


    #1.合并中心病例间 时空距离 小于合并阈值要求的 以病例为中心的时空立方体思路：
    #建立一个连接成分的二维数组，每一行代表一个连接成分，行内存该连接成分每个立方体在“超出阈值的时空立方体”数组中的下标。
    #遍历识别出的超出阈值的时空立方体，逐个与连接成分数组里的每个立方体的下标进行比较，看是否在阈值内。   发现一个在的，就意味着和这个连接成分相连通，就可以记录后直接搜索下一个连接成分了。
    #如果只与一个连接成分相连通，就归入；如果没有连通的连接成分，就新建一个连接成分。   如果出现与多个现有的连接成分是连通的，就合并这些连接成分，并归入；

    
    
    #Part2 统计在合并阈值范围内的 超出阈值的以病例为中心的时空立方体
    outbreaks_num = 0; #统计符合爆发条件的时空立方体个数，大于等于爆发个数
    outbteaks_csae_cube_Merged_Indexs=np.zeros((outbreaks_case_num,outbreaks_case_num+1),dtype=np.int32);  #(num,Index1,Index2,...)  
    #用来记录需要合并的 病例中心方格；每一行存的 是要合并到一起的 病例中心方格 在 outbteaks_csae_cube_Case_Indexs 中的下标
    #outbteaks_csae_cube_Merged_Indexs[i][0]为保留空间，用来记录该次聚集爆发点中要合并的病例中心方格的个数; 

    for i in range(outbreaks_case_num):
        #print(({i}));
        if 0 == i%100:
            print(f'{time.strftime("%Y/%m/%d %H:%M:%S",time.localtime())}   已经处理 {i} 个 超出阈值的 以病例为中心的时空立方体，标识出{outbreaks_num}个连接成分......');
        date=datetime.strptime(CaseData[Date_FieldName][outbteaks_csae_cube_Case_Indexs[i][1]],"%Y/%m/%d");
        inThreshold = np.zeros(outbteaks_csae_cube_numMax+1,dtype=np.int32);    #用来存 距离在合并阈值内的连接成分的下标，inThreshold[0]为个数。
        for j in range(outbreaks_num):
            #print(({i},{j}));
            for k in range(1,outbteaks_csae_cube_Merged_Indexs[j][0]+1):
                #print(({i},{j},{k}));
                if abs(CaseData[Lon_FieldName][outbteaks_csae_cube_Case_Indexs[i][1]] - CaseData[Lon_FieldName][outbteaks_csae_cube_Case_Indexs[outbteaks_csae_cube_Merged_Indexs[j][k]][1]]) < Lon_Threshold*Lon_interval:
                    if abs(CaseData[Lat_FieldName][outbteaks_csae_cube_Case_Indexs[i][1]] - CaseData[Lat_FieldName][outbteaks_csae_cube_Case_Indexs[outbteaks_csae_cube_Merged_Indexs[j][k]][1]]) < Lat_Threshold*Lat_interval:
                        date2 = datetime.strptime(CaseData[Date_FieldName][outbteaks_csae_cube_Case_Indexs[outbteaks_csae_cube_Merged_Indexs[j][k]][1]],"%Y/%m/%d");
                        if abs(date-date2).days <= Time_Threshold*Time_interval:
                            #与该连接成分中某病例在距离阈值范围内，意味着与该连接成分相连通
                            #print(f'({i},{j},{k}) ：满足阈值要求');
                            inThreshold[0] = inThreshold[0] +1;
                            inThreshold[inThreshold[0]] = j;
                            break;
        #print(f'第{i}个超出阈值的时空立方体, 与{inThreshold[0]}个已标识连接成分相连通');
        if inThreshold[0] == 0:
            #意味着没有与其连通的连接成分，所以新建一个连接成分
            outbteaks_csae_cube_Merged_Indexs[outbreaks_num][0] = outbteaks_csae_cube_Merged_Indexs[outbreaks_num][0]+ 1;
            outbteaks_csae_cube_Merged_Indexs[outbreaks_num][outbteaks_csae_cube_Merged_Indexs[outbreaks_num][0]] = i;
            outbreaks_num = outbreaks_num +1;
        elif inThreshold[0] == 1:
            #意味着只与一个连接成分相连通，所以归入
            outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0] = outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0] + 1;
            outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0]] = i;
        elif inThreshold[0] > 1:
            print(f'第{i}个超出阈值的时空立方体, 与{inThreshold[0]}个已标识连接成分相连通');
            #意味着与多个现有的连接成分是连通的，要合并这些连接成分，并归入
            for m in range(inThreshold[0],1):  #从倒数处理，往下标最小的合并（下标最小的会保留，不会被注销），这样可以保证注销已经合并的连接成分之后，还未合并的连接成分下标不变。   
                #合并 后面的 到 下标最小的
                outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0]+1:outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0] + outbteaks_csae_cube_Merged_Indexs[inThreshold[m]][0]+1] = outbteaks_csae_cube_Merged_Indexs[inThreshold[m]][1:outbteaks_csae_cube_Merged_Indexs[inThreshold[m]][0]+1];
                outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0] = outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0] + outbteaks_csae_cube_Merged_Indexs[inThreshold[m]][0];
                #注销已经归并了的
                for n in range(inThreshold[m],outbreaks_num):
                    outbteaks_csae_cube_Merged_Indexs[n] = outbteaks_csae_cube_Merged_Indexs[n+1];
                outbreaks_num = outbreaks_num -1;
            #合并完成了，归入合并后的
            outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0] = outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0] + 1;
            outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][outbteaks_csae_cube_Merged_Indexs[inThreshold[1]][0]] = i;





    #Part3 前面只是统计了该合并的立方体，下面开始真正合并病例
    outbteaks_Case_Indexs = np.zeros((outbreaks_num,case_in_cluster_numMax+1),dtype=np.int32);  #(num,Index1,Index2,...)  
    #用来记录每个聚集爆发事件里所有病例的Index，通过 ID_of_illness[Index]等即可获取到病例ID等信息。
    #outbteaks_Case_Indexs[i][0]为保留空间，用来记录该次聚集爆发事件中病例的个数; 

    #前面只是统计了该合并的立方体，下面开始真正合并病例
    print("合并处理中......")
    #case_index_all = [];
    #case_index_num = 0;
    for i in range(outbreaks_num):
        case_index_temp = [];
        for j in range(1,outbteaks_csae_cube_Merged_Indexs[i][0]+1):
            index_i = outbteaks_csae_cube_Merged_Indexs[i][j];
            case_index_temp = np.concatenate([case_index_temp,outbteaks_csae_cube_Case_Indexs[index_i][1:outbteaks_csae_cube_Case_Indexs[index_i][0]+1]]);  #concatenate 拼接并展开为1维
            _,index_unique_temp = np.unique(case_index_temp,return_index=True); 
            case_index_temp = case_index_temp[index_unique_temp];
        
            #case_index_all = np.concatenate([case_index_all,outbteaks_csae_cube_Case_Indexs[index_i][1:outbteaks_csae_cube_Case_Indexs[index_i][0]+1]]);
            #case_index_all = np.unique(case_index_all); 
        if case_index_temp.size > case_in_cluster_numMax:
            print(f'\n\n\nError: 第{i}个聚集暴发事件的病例数为{case_index_temp.size}, 超过设定的阈值case_in_cluster_numMax: {case_in_cluster_numMax},统计无法继续进行。\n请为函数传入一个更大的 case_in_cluster_numMa\n\n\n');
            return ;
            

        outbteaks_Case_Indexs[i][0] = case_index_temp.size;
        #case_index_num = case_index_num + case_index_temp.size;
        outbteaks_Case_Indexs[i][1:outbteaks_Case_Indexs[i][0]+1] = case_index_temp;
    #print("数值若相等，则表明统计值正确概率很大：",case_index_num,case_index_all.size);
    print(f'outbreaks_num: {outbreaks_num}');


    #2.合并有相同病例的 超出阈值的 以病例为中心的聚集爆发事件（当时空距离阈值大于等于1单位时，理论上第一种合并后，不会出现这种合并）
    for i in range(outbreaks_num-1,-1,-1):     #从倒数处理，往下标最小的合并（下标最小的会保留，不会被注销），这样可以保证注销已经合并的连接成分之后，还未合并的连接成分下标不变。   
        for j in range(i-1,-1,-1):
            case_index_temp = np.concatenate((outbteaks_Case_Indexs[i][1:outbteaks_Case_Indexs[i][0]+1],outbteaks_Case_Indexs[j][1:outbteaks_Case_Indexs[j][0]+1])); #concatenate 拼接并展开为1维
            case_index_temp = np.unique(case_index_temp); 
            if case_index_temp.size == (outbteaks_Case_Indexs[i][0] + outbteaks_Case_Indexs[j][0]):
                #两个聚集爆发事件的所有病例下标合并，去除重复值之后，个数与两个事件病例数之和相等，说明并没有重复值
                continue;
            elif case_index_temp.size < (outbteaks_Case_Indexs[i][0] + outbteaks_Case_Indexs[j][0]):
                #说明有重复值，则合并
                #print("合并");
                outbteaks_Case_Indexs[j][0]  = case_index_temp.size;
                outbteaks_Case_Indexs[j][1:outbteaks_Case_Indexs[j][0]+1] = case_index_temp[0:case_index_temp.size];
                #注销已经归并了的
                for n in range(i,outbreaks_num):
                    outbteaks_Case_Indexs[n] = outbteaks_Case_Indexs[n+1];
                outbreaks_num = outbreaks_num -1;

            else:
                print("Error合并重复值");


    print(f"完成 合并超出阈值的 以病例为中心的时空立方体，该部分耗时 {(time.perf_counter() - Process_Part_start_time)/60} 分钟......");
    print(f"总耗时 {(time.perf_counter() - Process_start_time)/60} 分钟......\n\n\n");
    Process_Part_start_time=time.perf_counter();




    #Part4 统计连接成分的指标，也就是最终的聚集爆发点
    outbreaks_num_stats = np.zeros(outbreaks_num_Threshold_big,dtype=np.int32);       #统计符合爆发条件的每个连接成分里病例个数的频率分布，只统计到20个
    outbreaks_final = np.zeros((outbreaks_num,14),dtype=np.double);  
    #(Time_year., Lon, Lat) , (start_Time_year. ,end_Time_year. ) , age_below_3_percentage,  num, 
    #(start_lat. , end_lat.) (start_lon. , end_lon.), (radius_lon. ,radius_lat.), radius

            

    #统计连接成分的指标，也就是最终的聚集爆发点
    print("\n\n\n统计聚集爆发点指标......");

    for i in range(outbreaks_num):
        if 0 == i%10:
            print(f'{time.strftime("%Y/%m/%d %H:%M:%S",time.localtime())}   ，已经处理 {i} 个聚集爆发事件......');
        
        #计算该聚集爆发点的时空重心
        date = datetime.strptime(CaseData[Date_FieldName][outbteaks_Case_Indexs[i][1]],"%Y/%m/%d")
        Time_center = 0.0;
        Lon_center=0.0;
        Lat_center=0.0;
        start_lat = 1e10;
        start_lon = 1e10;
        end_lat = -1e10;
        end_lon = -1e10;
        start_date_of_cluster = date;
        end_date_of_cluster = date;

        for j in range(1,outbteaks_Case_Indexs[i][0]+1):
            date2 = datetime.strptime(CaseData[Date_FieldName][outbteaks_Case_Indexs[i][j]],"%Y/%m/%d");
            if date2 < start_date_of_cluster:
                start_date_of_cluster = date2;
            if date2 > end_date_of_cluster:
                end_date_of_cluster = date2;
            if CaseData[Lon_FieldName][outbteaks_Case_Indexs[i][j]] < start_lon:
                start_lon = CaseData[Lon_FieldName][outbteaks_Case_Indexs[i][j]];
            if CaseData[Lon_FieldName][outbteaks_Case_Indexs[i][j]] > end_lon:
                end_lon = CaseData[Lon_FieldName][outbteaks_Case_Indexs[i][j]];
            if CaseData[Lat_FieldName][outbteaks_Case_Indexs[i][j]] < start_lat:
                start_lat = CaseData[Lat_FieldName][outbteaks_Case_Indexs[i][j]];
            if CaseData[Lat_FieldName][outbteaks_Case_Indexs[i][j]] > end_lat:
                end_lat = CaseData[Lat_FieldName][outbteaks_Case_Indexs[i][j]];
            Time_center = Time_center +(date2 - date).days;
            Lon_center += CaseData[Lon_FieldName][outbteaks_Case_Indexs[i][j]];
            Lat_center += CaseData[Lat_FieldName][outbteaks_Case_Indexs[i][j]];

        
        date = date + timedelta(days=Time_center/outbteaks_Case_Indexs[i][0]);
        outbreaks_final[i][0] = float(date.strftime("%Y.%m%d"));
        outbreaks_final[i][1] = Lon_center/outbteaks_Case_Indexs[i][0];
        outbreaks_final[i][2] = Lat_center/outbteaks_Case_Indexs[i][0];
        outbreaks_final[i][3] = float(start_date_of_cluster.strftime("%Y.%m%d"));
        outbreaks_final[i][4] = float(end_date_of_cluster.strftime("%Y.%m%d"));
        outbreaks_final[i][7] = start_lat;
        outbreaks_final[i][8] = end_lat;
        outbreaks_final[i][9] = start_lon;
        outbreaks_final[i][10] = end_lon;
        outbreaks_final[i][11] = (end_lon-start_lon)/2.0;
        outbreaks_final[i][12] = (end_lat-start_lat)/2.0;
        outbreaks_final[i][13] = np.sqrt(outbreaks_final[i][11]*outbreaks_final[i][12]);
        

        outbreaks_final[i][6] = outbteaks_Case_Indexs[i][0];

        if outbteaks_Case_Indexs[i][0] < outbreaks_num_Threshold_big:
            outbreaks_num_stats[outbteaks_Case_Indexs[i][0]] += 1;


        

    print(f"\n\n\n\n\n总计{outbreaks_num}个聚集爆发点(连接成分)");
    print("频率分布: ")
    for i in range(outbreaks_num_Threshold_big):
        print(f'{i}: {outbreaks_num_stats[i]};   累计：{outbreaks_num_stats[0:i+1].sum()}');

    print(f'病例数{outbreaks_num_Threshold_big}及以上的聚集爆发点：共{outbreaks_num - outbreaks_num_stats.sum()}个')
    for i in range(outbreaks_num):
        if outbteaks_Case_Indexs[i][0] >= outbreaks_num_Threshold_big:
            print(f'index: {i} ;   位置：({outbreaks_final[i][0]},{outbreaks_final[i][1]},{outbreaks_final[i][2]});   case_num: {outbteaks_Case_Indexs[i][0]} \n\n');


    print(f"完成病例时空密度统计, 总耗时 {(time.perf_counter() - Process_start_time)/60} 分钟......\n\n\n");

    return [outbreaks_num, outbreaks_final, outbteaks_Case_Indexs, outbreaks_num_stats, outbreaks_case_num, outbteaks_csae_cube_Case_Indexs, outbteaks_csae_cube_Merged_Indexs];



Time_interval = 7.0 ;  #3.0、7.0  单位为天
Lon_interval = 0.003 ; #单位为度
Lat_interval = 0.003 ;  #单位为度 
outbreaks_num_Threshold = 10;     #5、10  


Time_Threshold = 1;    # 1 倍 Time_interval
Lon_Threshold = 1;     # 1 倍 Lon_interval
Lat_Threshold = 1;     # 1 倍 Lat_interval
Date_FieldName="发病日期";
Lon_FieldName="lon"; 
Lat_FieldName="lat";   
case_in_cube_numMax = 1000;
case_in_cluster_numMax = 1000;
outbteaks_csae_cube_numMax = 10000;
outbreaks_num_Threshold_big = 20;


CaseData=pd.read_csv("CaseData_forTest.csv",sep=",");  #replace with your data

[outbreaks_num, outbreaks_final, outbteaks_Case_Indexs, outbreaks_num_stats, 
 outbreaks_case_num, outbteaks_csae_cube_Case_Indexs, outbteaks_csae_cube_Merged_Indexs] = Get_ST_Scan(CaseData,Time_interval, Lon_interval, Lat_interval, outbreaks_num_Threshold, 
                Time_Threshold, Lon_Threshold, Lat_Threshold, Date_FieldName, Lon_FieldName, Lat_FieldName, outbreaks_num_Threshold_big,
                case_in_cube_numMax, case_in_cluster_numMax, outbteaks_csae_cube_numMax);
