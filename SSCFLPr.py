# -*- coding: utf-8 -*-
## yfkong@henu.edu.cn, Apr.,2021
## An algorithm to solve problem of SSCFLP with service radius and covering percentage

import sys,os,random,time,copy,math,tempfile
#ArcGIS
has_arcpy=0
try:
    import arcpy
    has_arcpy=1
except:
    has_arcpy=0
#mip solver
mip_solvers=[] #MIP solvers supported 
mip_solver=''  #MIP solver, "cplex", "cbc" or ""
mip_file_path=tempfile.gettempdir()
os.chdir(mip_file_path)  #used in arcgis
try:
    import cplex
    #mip_solvers.append('cplex')
except: 
    pass
try:
    import pulp
    s=pulp.apis.GUROBI_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('gurobi')
    s=pulp.apis.cplex_api.CPLEX_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('cplex')
    s=pulp.apis.COIN_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('cbc')
except: 
    pass
if len(mip_solvers)>0: mip_solver=mip_solvers[0]

#constant
MAXNUMBER=1.0e+20
MINNUMBER=1.0e-10
#instance info
nodes=[]
nodes_std=[] #for homoDP only
weight_features=[] #for homoDP only
num_units=-1
nodedij=[]
nodedik=[]  #weighted cost from i to k, =nodedij*nodes[][3] 
nodendij=[] #network distance
node_neighbors=[]
facility_neighbors=[]
total_pop=0
avg_pop=0
total_supply=0
all_units_as_candadate_locations=1
facilityCandidate=[]
facilityCapacity=[]
facilityCost=[]
num_facilityCandidate=-1
num_districts=-1 # number of service areas/facilities
avg_dis_min=0.0
potential_facilities=[]
NearFacilityList=[]
nearCustomer=[]
nearCustomers=[]
geo_instance=1
pmp_I_eaquls_J=1

#parameters for districting
location_problem=1
max_num_facility=999
adaptive_number_of_facilities=1
fixed_cost_obj=1
spatial_contiguity=1 # 0 no, 1 yes, 2 yes with multi-parts
spatial_contiguity_minimum_percentage=5
pop_dis_coeff=10000.0 #used in the objective function
pop_deviation=0.00 #for pdp, 5%

#current solution
centersID=[]
node_groups=[]
district_info=[] #[[0,0,0.0] for x in range(num_districts)] # solution
objective_overload=0
obj_balance=MAXNUMBER
objective=MAXNUMBER
objective_fcost=MAXNUMBER
biobjective=MAXNUMBER
objective_supply=0.0

given_solution=0 #reserved
all_solutions=[]

#best solution in each start
best_solution =[] # node_groups[:]
best_centersID=[]
best_biobjective=MAXNUMBER
best_objective=MAXNUMBER
best_objective_overload = MAXNUMBER
best_objective_fcost = MAXNUMBER
#global best solution 
#best_centers_global=[]
best_solution_global=[]
best_centersID_global=[]
best_biobjective_global = MAXNUMBER
best_objective_global = MAXNUMBER
best_objective_fcost_global = MAXNUMBER
best_overload_global = MAXNUMBER

#search statistics
time_check=0
time_check_edge_unit=0
time_spp=0.0
time_update_centers=0.0
time_op=[0.0 for x in range(10)]
time_ruin_recreate=[0.0 for x in range(10)]
time_location=[0.0 for x in range(10)]
time_pmp_re_location=0.0
time_Whitaker=0.0
time_repair=0
count_op=[0.0 for x in range(10)]
check_count=0
improved=0
move_count=0

#search histry
region_pool = []
pool_index=[]

#local search
acceptanceRule="hc" #solver name
assignment_operators_selected=[0,1] #0 one-unit move, 1 two-unit move, 2 three-unit move
location_operators_selected=[0,1,2,3,4] #0 swap, 1 drop, 2 add, 3 add+drop, 4 me
ruin_oprators=[0,1,2,3,4] #ruin0, ruin1, 9 mip assign
multi_start_count=6 #population size for GA, ILS, VNS, LIS+VND
initial_solution_method=0 #0 construction, 1 LP
assign_method=0 #not used
assign_or_Location_search_method=0
large_facility_cost=0
maxloops=1000
max_loops_solution_not_improved=-1
SA_maxloops = 100 # maximum number of search loops for GA
SA_temperature=1.0
op_random = 1 # operators used sequentially (0) or randomly(1)
last_loop_obj=0.0
ruin_percentage=20
mainloop=0
mutation_rate=0.01 
cross_methods=-1
adj_pop_loops=5000
solution_similarity_limit=10.0

#mip modelling for inititail solution, spp and full model
is_spp_modeling=1 #0 no spp modelling, 1 modelling at the end, 2 modelling in case of local optimum
linear_relaxation=0
spp_loops=400
solver_time_limit=7200 #used for mip modeling
solver_mipgap=0.000000000001
solver_message=0
heuristic_time_limit=300
seed =random.randint(0,10000)
random.seed(seed)
locTabuLength=100  #nf*psize
locTabuList=[]
locTabuList2=[]

def update_locTabuList(locs):
    global locTabuList
    if locs not in locTabuList:
        locTabuList.append(locs)
    if len(locTabuList) >locTabuLength:
        del locTabuList[0]

def arcpy_print(s):
    if has_arcpy==1: 
        arcpy.AddMessage(s)
    else:
        print s

#record a region in current solution
def update_region_pool(rid):
    global pool_index
    global time_spp
    global region_pool
    t=time.time()
    if is_spp_modeling<=0: return
    if centersID[rid]<0: return
    #if spatial_contiguity==1 and check_continuality_feasibility(node_groups,rid)<1:
    #    return
    ulist=[x for x in  range(num_units) if node_groups[x]==rid]
    if ulist==[]:
        #print "empty area:",rid,node_groups
        return
    cost1=district_info[rid][2]
    cost2=sum(nodedik[x][rid] for x in ulist)
    if abs(cost1-cost2)>0.001: print rid,cost1,cost2
    obj=district_info[rid][2]+district_info[rid][4]*pop_dis_coeff
    idx=int(obj*100000)
    idx+=sum(x*x for x in ulist)
    if idx not in pool_index[rid]:
        pool_index[rid].append(idx)
        region_pool.append([ulist,district_info[rid][2],district_info[rid][1],district_info[rid][4],rid])
    time_spp+=time.time()-t
    return

#record all regions in current solution
def update_region_pool_all():
    if is_spp_modeling<=0:
        return
    for rid in range (num_districts):
        if centersID[rid]<0: continue
        update_region_pool(rid)

#check continuality of a solution (sol)
def check_solution_continuality_feasibility(sol):
    if spatial_contiguity==0: return 9
    feasible = spatial_contiguity
    for i in range (num_districts):
        if centersID[i]<0: continue
        if check_continuality_feasibility(sol,i) <= 0:
            feasible=0  #infeas.
            break
    return feasible

def check_current_solution_continuality_feasibility():
    feaslist=[]
    for i in range (num_districts):
        if centersID[i]<0: continue
        sta=check_continuality_feasibility(node_groups,i)
        feaslist.append(sta)
    return feaslist

#check continuality of a region (rid) in solution (sol)
def check_continuality_feasibility(sol,rid):
    global time_check
    global check_count
    if spatial_contiguity==0: return 9
    #if geo_instance==0: return -1
    u=facilityCandidate[rid]
    check_count+=1
    t=time.time()
    ulist1=[x for x in range(num_units) if sol[x]==rid and x!=u]
    ulist2=[u]
    #ulist2.append(ulist1.pop())
    for x in ulist2:
        for i in range(len(ulist1)):
            j=ulist1[i]
            if j in node_neighbors[x]:
                ulist2.append(j)
                ulist1[i]=-1
        ulist1=[x for x in ulist1 if x>=0]
    if len(ulist1)==0:          
        time_check+=time.time()-t
        return 1  #feasible
    if spatial_contiguity==2:
        ulist=[x for x in range(num_units) if sol[x]==rid]
        flist=frag_unit_minority(ulist)
        if flist==[]: 
            time_check+=time.time()-t
            return 2
    time_check+=time.time()-t
    return 0    #infeasible

#check continuality of a list of units
def check_ulist_continuality(ulist):
    if spatial_contiguity==0: return 1
    global time_check
    global check_count
    t=time.time()
    ulist1=ulist[:]
    ulist2=[]
    ulist2.append(ulist1.pop())
    check_count+=1
    for x in ulist2:
        for i in range(len(ulist1)):
            #if ulist1[i]==-1: continue
            if ulist1[i] in node_neighbors[x]:
                ulist2.append(ulist1[i])
                ulist1[i]=-1
        ulist1=[x for x in ulist1 if x>=0]         
    #ulist3=[x for x in ulist1 if x!=-1]
    if len(ulist1)==0:          
        time_check+=time.time()-t
        return 1  #feasible
    if spatial_contiguity==2:
        flist=frag_unit_minority(ulist)
        time_check+=time.time()-t
        if flist==[]: return 2
    return 0    #infeasible

#return a list of boundary units
def find_edge_units():
    if spatial_contiguity==0 and random.random()>0.5:
        return find_tail_units()
    ulist=[]
    for x in range(num_units):
        if geo_instance==1 and spatial_contiguity==1 and x in centersID:
            continue  #bug for benckmark instances
        k=node_groups[x]
        for y in node_neighbors[x]:
            if node_groups[y] != k:
                ulist.append(x)
                break
    random.shuffle(ulist)
    if objective_overload==0: 
        return ulist
    ulist=[[x, district_info[node_groups[x]][4]] for x in ulist]
    ulist.sort(key=lambda x:-x[1])
    ulist=[ x[0] for x in ulist]
    return ulist

#update region information of the current solution
def update_district_info():
    global objective_overload
    global objective
    global biobjective
    global objective_fcost
    global district_info
    global move_count
    global obj_balance
    global centersID
    global objective_supply
    global avg_dis_min
    for k in range(num_districts):
        district_info[k][0] = 0
        district_info[k][1] = 0.0
        district_info[k][2] = 0.0
        district_info[k][3] = 0.0
        district_info[k][4] = 0.0
    for k in range(num_districts):
        if centersID[k]<0 and k in node_groups:
            arcpy_print("debug: a facility not selected but used: " + str(k))
            centersID[k]=facilityCandidate[k]
    for k in range(num_districts):
        if centersID[k]<0:
            continue
        ulist=[x for x in range(num_units) if node_groups[x]==k]
        if len(ulist)==0:
            if location_problem==3: continue
            if adaptive_number_of_facilities==1:
                supply=sum(facilityCapacity[x] for x in range(num_districts) if centersID[x]>=0)
                centersID[k]=-1
                continue
        district_info[k][0] = len(ulist)
        district_info[k][1] = sum(nodes[x][3] for x in ulist)
        district_info[k][2] = sum(nodedik[x][k] for x in ulist)
        district_info[k][3] = facilityCapacity[k] 
        district_info[k][4] = max(0.0,district_info[k][1]-facilityCapacity[k]) # -district_info[k][3]
        if location_problem==3: district_info[k][4]=0 #pmp
        if location_problem==2: #pdp,edp
            bal=0.0
            dev=pop_deviation*total_pop/max_num_facility
            if district_info[k][1]>district_info[k][3]+dev: bal=district_info[k][1]-district_info[k][3]-dev
            if district_info[k][1]<district_info[k][3]-dev: bal=district_info[k][3]-district_info[k][1]-dev
            district_info[k][4]=bal
        #print centersID,node_groups
    bal=sum(x[4] for x in district_info)
    objective=sum([x[2] for x in district_info])
    objective_overload=bal
    #if objective/total_pop<avg_dis_min:
    #    avg_dis_min=objective/total_pop
    avg_dis_min=objective/total_pop
    biobjective=objective+objective_overload*avg_dis_min*pop_dis_coeff


    objective_supply=sum(facilityCapacity[x] for x in range(num_districts) if centersID[x] >=0)
    #biobjective=objective+objective_overload*avg_dis_min*1000000
    #biobjective=bal2*avg_dis_min*1000000
    if fixed_cost_obj==1:
        fcost=sum(facilityCost[x] for x in range(num_districts) if centersID[x] >=0)
        objective_fcost=fcost
        biobjective+=fcost
    move_count+=1

def find_frag_unit():
    if spatial_contiguity==0: return []
    global time_check_edge_unit
    t=time.time()    
    frag_units=[]
    if spatial_contiguity==2:
        for k in range(num_districts):
            if centersID[k]==-1: continue
            ulist=[x for x in range(num_units) if node_groups[x]==k]
            nflist=frag_unit_minority(ulist)
            frag_units+=nflist
    if spatial_contiguity!=2:
        for k in range(num_districts):
            if centersID[k]==-1: continue
            ulist2=[centersID[k]]
            ulist1=[x for x in range(num_units) if node_groups[x]==k and x!=centersID[k]]
            for x in ulist2:
                for i in range(len(ulist1)):
                    if ulist1[i]==-1: continue
                    if ulist1[i] in node_neighbors[x]:
                        ulist2.append(ulist1[i])
                        ulist1[i]=-1
                ulist1=[x for x in ulist1 if x>=0]
            frag_units+=ulist1

    random.shuffle(frag_units)
    time_check_edge_unit+=time.time()-t
    #print frag_units
    return frag_units    

def frag_unit_minority(ulist):
    final_list=[]
    ulist2=ulist[:1]
    ulist1=ulist[1:]
    total_area=sum(x[3] for x in nodes)
    while 1:
        for x in ulist2:
            for i in range(len(ulist1)):
                if ulist1[i]==-1: continue
                if ulist1[i] in node_neighbors[x]:
                    ulist2.append(ulist1[i])
                    ulist1[i]=-1
            ulist1=[x for x in ulist1 if x>=0]
        final_list.append([len(ulist2),ulist2[:]])
        if len(ulist1)<=1:
           if len(ulist1)==1:
                final_list.append([1,ulist1[:]])
           break
        u=ulist1[0]
        ulist2=[u]
        del ulist1[0]
    if len(final_list)==1: return []
    final_list.sort(key=lambda x:x[0])
    #del final_list[-1]
    flist=[]
    n=total_area*spatial_contiguity_minimum_percentage/max_num_facility/100
    for x in final_list: 
        area=sum(nodes[i][3] for i in x[1])
        if area>n:continue
        flist+=x[1]
    #print [len(flist),len(ulist)],
    return flist

def district_with_multiparts(ulist):
    final_list=[]
    ulist2=ulist[:1]
    ulist1=ulist[1:]
    while 1:
        for x in ulist2:
            for i in range(len(ulist1)):
                if ulist1[i]==-1: continue
                if ulist1[i] in node_neighbors[x]:
                    ulist2.append(ulist1[i])
                    ulist1[i]=-1
            ulist1=[x for x in ulist1 if x>=0]
        final_list.append([len(ulist2),ulist2[:]])
        if len(ulist1)<=1:
           if len(ulist1)==1:
                final_list.append([1,ulist1[:]])
           break
        u=ulist1[0]
        ulist2=[u]
        del ulist1[0]
    final_list.sort(key=lambda x:-x[0])
    #del final_list[-1]
    flist=[x[0] for x in final_list]
    return flist

def repair_fragmented_solution_large_instance():
    global node_groups
    global centersID
    global time_repair
    if spatial_contiguity==0: return
    t=time.time()
    frag_units=find_frag_unit()
    #print "frag_units",len(frag_units),
    if len(frag_units)==0: return
    #sol=node_groups[:]
    for x in frag_units:
        node_groups[x]=-1
    update_district_info()
    #print len(frag_units),
    while len(frag_units)>0:
        cands=[]
        for x in frag_units:
            for y in node_neighbors[x]:
                k=node_groups[y]
                if k<0: continue
                cands.append([x,k,nodedik[x][k]])
        if cands==[]: break
        cands.sort(key=lambda x:x[2])
        n=len(cands)/10
        if n<1: n=1
        for x in cands[:n]:
            nid=x[0]
            if node_groups[nid]>=0: continue
            node_groups[nid]=x[1]
            if nid in frag_units: frag_units.remove(x[0])
    for x in frag_units:
        for k in NearFacilityList[x]:
            if centersID[k]>=0:
                 node_groups[x]=k
                 break
    #print len(frag_units)
    #print int(objective),
    update_district_info()
    #print int(objective),
    time_repair+=time.time()-t
    repair_fragmented_solution()
    #print objective,

#repair the fragmented solution
def repair_fragmented_solution():
    if spatial_contiguity==0: return
    #if num_units/max_num_facility>100:
    #    repair_fragmented_solution_large_instance()
    #    return
    global node_groups
    global centersID
    global time_repair
    t=time.time()
    for k in range(num_districts):
        if centersID[k]<0: continue
        u=nearCustomer[k]
        node_groups[u]=k
    update_district_info()
    frag_units=find_frag_unit()
    #print "frag_units",frag_units,
    if len(frag_units)==0: return
    sol=node_groups[:]
    for x in frag_units:
        node_groups[x]=-1
    # if location_problem>=3:
        # for k in range(num_districts):
            # if centersID[k]<0: continue
            # if node_groups[k]>=0: continue
            # c,ulist=update_center(k)
            # centersID[k]=-1
            # centersID[c]=c
            # for x in ulist: node_groups[x]=c
    update_district_info()
    #print len(frag_units),
    while len(frag_units)>0:
        newk=-1
        nid=-1
        cost=MAXNUMBER
        for x in frag_units:
            for y in node_neighbors[x]:
                k=node_groups[y]
                if k>=0:
                    gap=max(district_info[k][1]+ nodes[x][3] - facilityCapacity[k],0)
                    if location_problem==3: gap=0
                    cost2=gap*avg_dis_min*pop_dis_coeff + nodedik[x][k]
                    if cost2<cost:
                        nid=x
                        newk=k
                        cost=cost2
        if newk>=0:
            node_groups[nid]=newk
            update_district_info()
            frag_units.remove(nid)
        else:
            break
    #print "frag_units", frag_units
    for x in frag_units:
        node_groups[x]=sol[x]
    #print len(frag_units)
    #print int(objective),
    update_district_info()
    #print check_current_solution_continuality_feasibility()
    time_repair+=time.time()-t
    #print int(objective),

def select_region(seed):
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    n=100*nf/num_units  #areas with 100 units
    if nf<=5: n=nf
    #if nf>=7: n=random.randint(nf/2+1,nf)
    if nf>=6 and nf<=11: n=random.randint(nf/2+1,9)
    if nf>=12 and nf<=16: 
        n=random.randint(7,10)
    if nf>=17: 
        n=random.randint(7,10)
    if n*num_units/nf<80: 
        n=min(10,80*nf/num_units)
    if location_problem==3: n=min(128,num_districts/10)
    #clist=[]
    #u=random.randint(0,num_units-1)
    #if r>=0: 
    #    u=nearCustomer[r]
    #    clist=[r]
    #for k in NearFacilityList[u]:
    #    if centersID[k]<0: continue
    #    if k not in clist: clist.append(k)
    #    if len(clist)==n: break
    #return clist
    #if objective_overload>0: ???
    u=seed
    if u<0: u=random.randint(0,num_units-1)
    r=node_groups[u]
    if location_problem==0 and objective_overload>0: #SAP
        rlist=[k for k in range(num_districts) if district_info[k][4]>0]
        random.shuffle(rlist)
        r=rlist[0]
        u=nearCustomer[r]
        return NearFacilityList[u][:n]

    clist=[r]
    if random.random()>-0.5:
        for k in NearFacilityList[u]:
            if centersID[k]<0: continue
            if k==r: continue
            clist.append(k)
            if len(clist)==n: break
        #clist.sort()
        return clist

    for i in facility_neighbors[r]:
        if centersID[i]<0: continue
        clist.append(i)
        if len(clist)>=n: break
    #clist.sort()
    return clist

def select_region_cflpr(seed):
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    if nf<=10: n=random.randint(nf/2+1,nf)
    if nf>=11 and nf<=20:
         n=random.randint(9,11)
    if nf>=21: 
        n=random.randint(10,13)
    u=random.randint(0,num_units-1)
    r=node_groups[u]
    clist=[r]
    for k in NearFacilityList[u]:
        if centersID[k]<0: continue
        if k==r: continue
        clist.append(k)
        if len(clist)==n: break
    return clist


#update the best and the global best solution
#if the current solution is better than the best
def update_best_solution():
    global best_solution
    global best_centersID
    global best_biobjective
    global best_objective
    global best_objective_fcost
    global best_overload
    global best_objective_overload
    global best_centersID_global
    global best_solution_global
    global best_biobjective_global
    global best_objective_global
    global best_objective_fcost_global
    global best_overload_global    
    global improved_loop
    global improved
    global avg_dis_min
    improve =0
    if location_problem==1 and adaptive_number_of_facilities==0:
        nf=sum(1 for x in centersID if x>=0)
        if nf!=max_num_facility: return 0
    #if spatial_contiguity==1 and check_solution_continuality_feasibility(node_groups)==0:
    #    ##noprint "check_solution_continuality_feasibility!!!"
    #    return improve
    biobj=biobjective
    biobj_best=best_biobjective
    if biobj<=biobj_best:
        best_biobjective=biobj
        best_objective = objective
        best_objective_fcost=objective_fcost
        best_objective_overload = objective_overload
        best_solution = node_groups[:]
        best_centersID=centersID[:]
        improved_loop=mainloop
        improve =1
        improved+=1
    if biobj<best_biobjective_global:
        best_biobjective_global=biobj
        best_centersID_global=centersID[:]
        best_overload_global = objective_overload
        best_solution_global =node_groups[:]
        best_objective_global = objective
        best_objective_fcost_global=objective_fcost
        avg_dis_min=biobj/total_pop
    return improve

#return the neighor regions of unit nid
def neighbor_regions(nid):
    rid=node_groups[nid]
    if spatial_contiguity>=1:
        rlist2=[node_groups[x] for x in node_neighbors[nid] if node_groups[x]!=rid]
        rlist=list(set(rlist2))
        if len(rlist2)>1: random.shuffle(rlist2)
        return rlist
    if spatial_contiguity>=0 and random.random()>0.5: #testing ??? 
        knn=4
        #if knn< num_units/100: knn=num_units/100
        rlist=[]
        for k in NearFacilityList[nid]:
            if k==rid: continue
            if centersID[k]<0: continue
            rlist.append(k)
            if len(rlist)>=knn: break 
        return rlist
    rlist2=[node_groups[x] for x in node_neighbors[nid] if node_groups[x]!=rid]
    rlist=list(set(rlist2))
    if len(rlist2)>1: random.shuffle(rlist2)
    return rlist

def create_facility_neighbors():
    return 
    global facility_neighbors
    mindij=[[MAXNUMBER for x in range(num_districts)] for y in range(num_districts)]
    for i in range(num_districts):
        for j in range(num_districts):
            if j<=i: continue
            dlist=[nodedij[x][i]-nodedij[x][j] for x in range(num_units)]
            d=sum(x*x for x in dlist)
            mindij[i][j]=d
            mindij[j][i]=d
    facility_neighbors = [[]for x in range(num_districts)]
    for i in range(num_districts):
        dlist=[[x, mindij[i][x]] for x in range(num_districts) ]
        dlist.sort(key=lambda x:x[1])
        nghrs=[x[0] for x in dlist]
        facility_neighbors[i]=nghrs[:]

def create_node_neighbors():
    global node_neighbors
    #rlist=[x for x in range(num_districts)]
    mindij=[[MAXNUMBER for x in range(num_units)] for y in range(num_units)]
    for i in range(num_units):
        for j in range(num_units):
            if j<=i: continue
            dlist=[nodedij[i][x]-nodedij[j][x] for x in range(num_districts)]
            d=sum(x*x for x in dlist)
            mindij[i][j]=d
            mindij[j][i]=d
    node_neighbors = [[]for x in range(num_units)]
    for i in range(num_units):
        dlist=[[x, mindij[i][x]] for x in range(num_units) ]
        dlist.sort(key=lambda x:x[1])
        nn=8
        if nn>num_units: nn=num_units
        nghrs=[dlist[x][0] for x in range(nn)]
        random.shuffle(nghrs) #debug
        node_neighbors[i]=nghrs[:]

#read instance file, f1:unit info, f2: connectivity info
def readfile(f1):
    global num_units
    global total_pop
    global total_supply
    global nodes
    global node_neighbors
    global nodedij
    global nodedik
    global centersID
    global facilityCandidate
    global facilityCapacity
    global facilityCost
    global num_facilityCandidate
    global num_districts
    global district_info
    global avg_dis_min
    global potential_facilities
    node =[0,0,0,0,0,0,0,0,0,0]
    #school_nodes = []
    nodes = []
    #nodes.append(node)
    ##noprint "reading nodes ...",
    f = open(f1)
    line = f.readline()  #OID    pop    PointX    PointY    fcadidature    fcost    fcapacity
    line = f.readline()
    nodeidx=0
    while line:
        line=line[0:-1]
        #print line
        items = line.split(',')
        if len(items)<=2:
            items = line.split('\t')
        unit=[nodeidx, float(items[2]), float(items[3]), int(items[1]),int(items[0]),int(items[6]),int(items[4]),float(items[5])]
        nodes.append(unit)
        nodeidx+=1
        #nodes.append([int(items[1]), float(items[8]), float(items[9]), int(items[5]), int(items[6]), int(items[7]),int(items[12]),int(items[13])])
        line = f.readline()
    f.close()
    num_units=len(nodes)
    total_pop=sum(x[3] for x in nodes)
    ##noprint num_units,"units"
    ##noprint "reading connectivity ...",
    connectivity=[[0 for x in range(len(nodes))] for y in range(len(nodes))]
    ###id1,id2#####
    '''
    f = open(f2)
    line = f.readline()
    line = f.readline()
    links=0
    while line:
        items = line.split(',')
        if len(items)<=2:
            items = line.split('\t')
        if int (items[1]) != int (items[2]):
            id1=int (items[1])
            id2=int (items[2])
            idx1=-1
            idx2=-1
            for i in range(num_units):
                if nodes[i][4]==id1:
                    idx1=i
                if nodes[i][4]==id2:
                    idx2=i
                if idx1>=0 and idx2>0:
                    break
            connectivity[idx1][idx2]=1
            connectivity[idx2][idx1]=1
            links+=1
        line = f.readline()
    f.close()'''
    ##noprint links,"links"
    num_units=len(nodes)
    facilityCandidate=[]
    facilityCapacity=[]
    facilityCost=[]
    centersID=[]
    ##noprint "all data are read! "
    for i in range(num_units):
        if nodes[i][5]>0 or all_units_as_candadate_locations==1:
            facilityCandidate.append(i)
            facilityCapacity.append(nodes[i][5])
            facilityCost.append(nodes[i][7])
            centersID.append(i)
    num_facilityCandidate=len(facilityCandidate)
    num_districts=len(facilityCandidate)
    #facilityCandidate.sort()
    total_supply=sum(facilityCapacity)
    centersID=facilityCandidate[:]
    nodedij=[[MAXNUMBER for x in range(num_districts)] for y in range(num_units)]
    max_dij=0.0
    for i in range(num_units):
        for k in range(num_districts):
            j=facilityCandidate[k]
            d2=pow(nodes[i][1]-nodes[j][1],2)
            d2+=pow(nodes[i][2]-nodes[j][2],2)
            d=pow(d2,0.5)/1000
            nodedij[i][k]=d
            if d>max_dij:
                max_dij=d

    node_neighbors = [[]for x in range(len(nodes))]
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if j<=i:
                continue
            if connectivity[i][j]==1:
                node_neighbors[i].append(j)
                node_neighbors[j].append(i)

    district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
    dis=0.0
    for i in range(num_units):
        d=min(nodedij[i])
        dis+=d*nodes[i][3]
    avg_dis_min=dis/total_pop
    #weighted cost from i to k
    
    nodedik=[[nodedij[y][x]*nodes[y][3] for x in range(num_districts)] for y in range(num_units)]
    find_NearFacilityList(num_districts)
    print "find_near_customer()..."
    find_near_customer()
    #find_nearFacilityFacility()
    print "create_facility_neighbors()..."
    create_facility_neighbors()
    potential_facilities=[x for x in range(num_districts)]
    s="M N: "+str(num_districts)+" "+str(num_units)
    arcpy_print(s)
    s="Total demand & supply: "+str(total_pop)+" "+str(total_supply)
    arcpy_print(s)

def find_nearFacilityFacility():
    global nearFacilityFacility
    nearFacilityFacility=[[] for x in range(num_districts)]
    dkk=[sum(nodedik[x][k]*nodedik[x][k] for x in range(num_units)) for k in range(num_districts)]
    #dkk.sort(key=lambda x:x[1])
    #dlist=[x[0] for x in dkk]
    for k in range(num_districts):
        d=dkk[k]
        dk=[[x,dkk[x]-d] for x in range(num_districts)]
        dk.sort(key=lambda x:x[1])
        del dk[0]
        nearFacilityFacility[k]=[x[0] for x in dk]
def find_near_customer():
    global nearCustomer
    global nearCustomers
    if location_problem>=2 and pmp_I_eaquls_J==1: 
        nearCustomers=NearFacilityList
        nearCustomer=[x for x in range(num_districts)]
        return
    nearCustomer=[-1 for x in range(num_districts)]
    nearCustomers=[[] for x in range(num_districts)]
    dis=[]
    for k in range(num_districts):
        dis=[ [x,nodedij[x][k]] for x in range(num_units)]
        dis.sort(key=lambda x: x[1])
        nearCustomer[k]=dis[0][0]
        nearCustomers[k]=[x[0] for x in dis]
       
def initialize_instance():
    global num_units
    global num_districts
    global num_facilityCandidate
    global centersID
    global node_groups
    global facilityCost
    global facilityCandidate
    global facilityCapacity
    global nodedik
    global avg_pop
    global total_pop
    global avg_dis_min
    global total_supply
    global fixed_cost_obj
    global max_num_facility
    global adaptive_number_of_facilities
    #solution obj 
    global district_info
    global objective_overload
    global objective
    global biobjective
    global all_solutions
    #best solution 
    global best_solution
    global best_centersID
    global best_biobjective
    global best_objective
    global best_objective_overload
    #global best solution 
    global best_solution_global
    global best_centersID_global
    global best_biobjective_global
    global best_objective_global
    global best_overload_global
    global potential_facilities
    global max_exclusion_list
    num_districts=len(facilityCandidate)
    #num_units=len(nodes)
    total_pop=sum(x[3] for x in nodes)
    #print total_pop,nodes[:10]
	#sum(nodes[x][3] for x in range(num_units))
    node_groups=[-1 for x in range(num_units)]
    if location_problem==0:
        fixed_cost_obj=0
        max_num_facility=num_districts
    if fixed_cost_obj==0:
        facilityCost=[0 for x in range(num_districts)]
    if location_problem==1 and max_num_facility<1:
        max_num_facility=num_districts
        adaptive_number_of_facilities=1
    if location_problem==2:
        if all_units_as_candadate_locations==1:
            facilityCandidate=[x for x in range(num_districts)]
            facilityCost=[0.0 for x in range(num_districts)]
            popa=total_pop*1.0/max_num_facility
            facilityCapacity=[popa for x in range(num_districts)]
        if all_units_as_candadate_locations==0:
            facilityCost=[0.0 for x in range(num_districts)]
            popa=total_pop*1.0/max_num_facility
            facilityCapacity=[popa for x in range(num_districts)]
    if location_problem==3: #pmp
        #num_districts=num_units
        #facilityCandidate=[x for x in range(num_districts)]
        facilityCost=[0.0 for x in range(num_districts)]
        facilityCapacity=[total_pop for x in range(num_districts)]
    centersID=facilityCandidate[:]
    num_facilityCandidate=len(facilityCandidate)
    district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
    total_supply=sum(facilityCapacity)
    #arcpy_print("total demand: "+str(total_pop))
    #arcpy_print("total supply: "+str(total_supply))
    #arcpy_print("avg. distance to nearest facility: "+str(avg_dis_min))

    objective_overload=MAXNUMBER
    obj_balance=MAXNUMBER
    objective=MAXNUMBER
    biobjective=MAXNUMBER
    all_solutions=[]

    #best solution in each start
    best_solution =[] # node_groups[:]
    best_centersID=[]
    best_biobjective=MAXNUMBER
    best_objective=MAXNUMBER
    best_objective_overload = MAXNUMBER

    #global best solution 
    best_solution_global=[]
    best_centersID_global=[]
    best_biobjective_global = MAXNUMBER
    best_objective_global = MAXNUMBER
    best_overload_global = MAXNUMBER
    #if geo_instance==1:
    #    nodedik=[[nodedij[y][facilityCandidate[x]]*nodes[y][3] for x in range(num_districts)] for y in range(num_units)]
    avg_dis_min =sum(nodedik[x][0] for x in range(num_units))/total_pop
    if spatial_contiguity>=1:
        find_near_customer()
    find_NearFacilityList(num_districts)
    if linear_relaxation==1:
        lplocs,sol=location_model_linear_relexation(max_num_facility,0,heuristic_time_limit,0.0001)	
        potential_facilities=[x for x in range(num_districts) if lplocs[x]>0.0001]
        print "Potential facilities by Linear Relax",potential_facilities    
    max_exclusion_list=[0.0 for x in range(num_districts)]

    
#read network distance
def readdistance(dfile):
    global nodedij
    nodedij=[[MAXNUMBER for x in range(len(nodes))] for y in range(len(nodes))]
    ##noprint "reading distances ...",
    try:
        f = open(dfile)
        line = f.readline()
        line = f.readline()
        readsuccess=1
        while line:
            items = line.split(',')
            if len(items)<=2:
                items = line.split('\t')
            if int (items[1]) != int (items[2]):
                id1=int (items[1])
                id2=int (items[2])
                idx1=-1
                idx2=-1
                for i in range(num_units):
                    if nodes[i][4]==id1:
                        idx1=i
                    if nodes[i][4]==id2:
                        idx2=i
                    if idx1>=0 and idx2>0:
                        break
                if idx1>=0 and idx2>=0:
                    nodedij[idx1][idx2]=float(items[3])
            line = f.readline()
        f.close()
        find_NearFacilityList(num_districts)
        return 1
    except:
        arcpy_print("Cannot read the distance data file!!!")
        return 0

def find_NearFacilityList(nnn):
    global NearFacilityList
    if len(NearFacilityList)>0: return
    NearFacilityList=[]
    n=nnn#num_districts
    if n>num_districts: n=num_districts
    dis=0.0
    print "find_NearFacilityList()",
    for i in range(num_units):
        if i%100==0: print ".",
        fdlist=[ [x,nodedik[i][x]] for x in range(num_districts)]
        fdlist.sort(key=lambda x:x[1])
        flist=[x[0] for x in fdlist[0:n]]
        NearFacilityList.append(flist[:])
    print
    if geo_instance==0:
        return
            
def get_cand_locations(dlist,ulist):
    clist=[]
    if len(max_exclusion_list)<=0: 
        clist=[x for x in range(num_districts) if centersID[x]<0 and nearCustomer[x] in ulist]
        random.shuffle(clist)
    else:
        clist=[[x,max_exclusion_list[x]] for x in range(num_districts) if centersID[x]<0 and nearCustomer[x] in ulist and x in potential_facilities]
        clist.sort(key=lambda x:x[1])
        clist=[x[0] for x in clist]
    nc=len(clist)
    nf=len(dlist)
    nu=len(ulist)
    t=time.time()
    mnsize=min(3000,num_units*num_districts/10)
    ns=min(nc,mnsize/nu-nf)
    #ns=nc
    if ns<nf: ns=min(nc,nf)
    if ns==nc: return [clist]
    #cloclist=clist[:ns/2]
    #while len(cloclist)<ns:
    #    r=random.random()
    #    idx=int(r**1.5*0.999*(nc-ns/2))
    #    k=clist[ns/2+idx]
    #    if k not in cloclist: cloclist.append(k)
    #cloclist.append(clist[:ns])
    cloclist=[]
    clist1=clist[:ns/2]
    clist2=clist[ns/2:]
    for i in range(5):
        if initial_solution_method!=9:
            random.shuffle(clist)
            cloclist.append(clist[:ns])
        else:
            random.shuffle(clist2)
            cloclist.append(clist1+clist2[:ns/2])
    #print nf,nu,mnsize,nc+nf,ns+nf,
    return cloclist

#capacitated location set covering 
def cflpr_model(max_radius,cover_pecentage,numf,time_limit,mipgap): 
    global centersID
    global node_groups
    global all_solutions
    prob = pulp.LpProblem("pcp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    xvariables={}
    costs={}
    ulist=range(num_units)
    for i in ulist:
        for j in centers:
            #if nodedij[i][j]>max_radius: continue
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        #costs["y_" +str(i)]=facilityCost[i]
    zvariable=pulp.LpVariable("z", 0, None, pulp.LpContinuous)
    fvariables={}
    if spatial_contiguity==1:
        for i in range(num_units):
            for j in node_neighbors[i]:
                for k in range(num_districts):
                    fvariables["f_" +str(i)+ "_"+str(j)+ "_"+ str(k)]=pulp.LpVariable("f_" +str(i)+ "_"+ str(j)+ "_"+ str(k), 0, None, pulp.LpInteger)

    obj=""
    for x in yvariables:
        obj +=100000000* yvariables[x]
    for i in ulist:
         for j in centers:
            obj+=costs["x_" +str(i)+ "_"+ str(j)]*xvariables["x_" +str(i)+ "_"+ str(j)]
    #obj+=100000000*zvariable
    prob += obj

    #cons 2
    for i in ulist:
        s=""
        for j in centers:
            #if cover_pecentage>=100 and nodedij[i][j]>max_radius: continue
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        prob += s==1

    for k in centers:
        s=""
        for i in ulist:
            ###if nodedij[i][k]>max_radius: continue
            s+= nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(k) ]
        s-= facilityCapacity[k]*yvariables["y_" +str(k)]
        prob += s <= 0
    
    if cover_pecentage<100:
        s=""
        for i in ulist:
            for j in centers:
                if nodedij[i][j]>max_radius: continue
                s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
        prob += s >=total_pop*cover_pecentage/100

    if numf>0:
        s=""
        for j in centers:
            s+=yvariables["y_" +str(j)]
        prob += s==numf
    #prob.writeLP("_pcp.lp")
    initvalues=cflp_mst()
    for x,v in initvalues:
       #items=x.split("_")
       if x.find("y")==0:  yvariables[x].setInitialValue(v)       
       if x.find("x")==0:  xvariables[x].setInitialValue(v)

    if spatial_contiguity==1:
        for i in range(num_units):
            for j in node_neighbors[i]:
                for k in range(num_districts):
                    s=fvariables["f_" +str(i)+ "_"+str(j)+ "_"+ str(k)]
                    s-= num_units * xvariables["x_" +str(i)+ "_"+ str(k)]
                    prob += s <= 0
        for i in range(num_units):
            for j in node_neighbors[i]:
                for k in range(num_districts):
                    s=fvariables["f_" +str(i)+ "_"+str(j)+ "_"+ str(k)]
                    s-= num_units* xvariables["x_" +str(j)+ "_"+ str(k)]
                    prob += s <= 0
        for i in range(num_units):
            for k in range(num_districts):
                if facilityCandidate[k]!=i:
                    s=""
                    for j in node_neighbors[i]:
                        s+=fvariables["f_" +str(i)+ "_"+str(j)+ "_"+ str(k)] - fvariables["f_" +str(j)+ "_"+str(i)+ "_"+ str(k)]
                    s-= xvariables["x_" +str(i)+ "_"+ str(k)]
                    prob += s >= 0

    #warmStart=True,
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)

    if prob.status<=0:
        print "model unsolved... or infeasible"
        return []
    centersID=[-1 for x in range(num_districts)]
    if len(node_groups)<1: node_groups=[-1 for x in range(num_units)]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            if v.name=="z": continue
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
            if items[0]=='x':
                node_groups[i]=int(items[2])
    update_district_info()
    coverdemand=0
    for i in range(num_units):
        k=node_groups[i]
        if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
    print "model solution:",biobjective,objective_fcost,objective,objective_overload,int(coverdemand*10000/total_pop)/100.0
    all_solutions.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost])
    all_solutions.sort(key=lambda x:x[0])
    return 1

def lscp_mip_model_init(max_radius,cover_pecentage,time_limit,mipgap): 
    global centersID
    global node_groups
    prob = pulp.LpProblem("pcp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    sampling=0
    if total_supply>15*total_pop:
        sampling=1
        min_nf=total_pop*num_districts/total_supply
        centers=k_medoids_cflpr_sampling(min_nf*15)
    if total_supply>15000*total_pop:
        centers=[]
        ilist=[]
        supply=0
        try_count=0
        while supply<15*total_pop:
            if try_count>num_districts: break
            try_count+=1
            k=random.randint(0,num_districts-1)
            if len(ilist)>0:
                mind=min(nodedij[x][k] for x in ilist)
                if mind<max_radius*0.5: continue
            if k not in centers:
                centers.append(k)
                ilist.append(nearCustomer[k])
                supply+=facilityCapacity[k]
                try_count=0
            #print [try_count,len(centers),supply],
    ulist=range(num_units)
    xvariables={}
    costs={}
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        costs["y_" +str(i)]=100000000#facilityCost[i]
    zvariables={}
    for i in ulist:
        zvariables["z_" +str(i)]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpContinuous)

    obj=""
    for x in yvariables:
        obj +=yvariables[x]
    prob += obj

    min_nf=total_pop*num_districts/total_supply
    #con 1
    #s=""
    #for k in centers:
    #    s+=yvariables["y_" +str(k)]
    #prob +=s >= min_nf

    s=""
    for k in centers:
        s+=facilityCapacity[k]*yvariables["y_" +str(k)]
    prob +=s >= total_pop

    #cons 2
    for i in ulist:
        s=""
        for j in centers:
            if nodedij[i][j]>max_radius: continue
            s+=yvariables["y_" + str(j)]
        prob +=s - zvariables["z_" +str(i)] >= 0

#    for i in ulist:
#        for j in centers:
#            if nodedij[i][j]>max_radius: continue
#            s=yvariables["y_" + str(j)]+ zvariables["z_" +str(i)]
#            prob +=s <= 1

    s=""
    for i in ulist:
        r=1
        if sampling==0:
            r=(49.5+random.random())/50
        s+=nodes[i][3]*r*zvariables["z_" +str(i)]
    prob += s>=total_pop*cover_pecentage/100
    #maxuc=total_pop-total_pop*cover_pecentage/100
    #prob += s<=maxuc
    prob.writeLP("_lscp.lp")
    #initvalues=pmp_mst(dlist,ulist)
    #for x,v in initvalues:
    #    if x.find("x")==0:  xvariables[x].setInitialValue(v)
    #    if x.find("y")==0:  yvariables[x].setInitialValue(v)
    #warmStart=True,
    #solver_message=1
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)

    if prob.status<=0:
        print "model unsolved..."
        return prob.status<0
    centersID=[-1 for x in range(num_districts)]
    if len(node_groups)<1: node_groups=[-1 for x in range(num_units)]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k] >=0:
                node_groups[i]=k
                break
    update_district_info()
    return 1

def clscp_mip_model_init(max_radius,cover_pecentage,time_limit,mipgap): 
    global centersID
    global node_groups
    prob = pulp.LpProblem("pcp",pulp.LpMinimize)
    centers=range(num_districts)#facilityCandidate
    sampling=0
    if total_supply>15*total_pop:
        sampling=1
        min_nf=total_pop*num_districts/total_supply
        centers=k_medoids_cflpr_sampling(min_nf*15)
        #print len(centers),centers
        #centers=[]
        #ilist=[]
        #supply=0
        #try_count=0
        #while supply<10*total_pop:
        #    if try_count>num_districts: break
        #    try_count+=1
        #    k=random.randint(0,num_districts-1)
        #    if len(ilist)>0:
        #        mind=min(nodedij[x][k] for x in ilist)
        #        if mind<max_radius*0.5: continue
        #    if k not in centers:
        #        centers.append(k)
        #        ilist.append(nearCustomer[k])
        #        supply+=facilityCapacity[k]
        #        try_count=0
        #    #print [try_count,len(centers),supply],
    covered_units=[]
    unit_covered=[0 for x in range(num_units)]
    for k in range(num_districts):
        covered=[]
        if k in centers:
            supply=facilityCapacity[k]*120/100  #pie*r*r / 2.6*r*r = pie/2.6=
            for i in nearCustomers[k]:
                if nodedij[i][k]>max_radius: break
                if supply<nodes[i][3]: break
                supply-=nodes[i][3]
                covered.append(i)
                unit_covered[i]=1
        covered_units.append(covered[:])

    #maxf=total_supply/total_pop
    #for i in range(num_units):  
    #    if unit_covered<maxf: continue
    #    for k in NearFacilityList[i]:
    #        if k in centers:
    #            covered_units[k].append(i)
    #        if len(covered_units[k])>=maxf: break

    ulist=range(num_units)
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
    zvariables={}
    for i in ulist:
        zvariables["z_" +str(i)]=pulp.LpVariable("z_" +str(i), 0, 1, pulp.LpContinuous)
    obj=""
    for x in yvariables:
        obj +=yvariables[x]
    prob += obj
    """
    min_nf=total_pop*num_districts/total_supply
    s=""
    for k in centers:
        s+=yvariables["y_" +str(k)]
    if adaptive_number_of_facilities==0 and max_num_facility>0: 
        prob +=s == max_num_facility
    else:
        prob +=s >= min_nf

    for i in ulist:
        s=""
        for j in centers:
            if i not in covered_units[j]: continue
            s+=yvariables["y_" + str(j)]
        prob += s>= 1
    """
    s=""
    for k in centers:
        s+=facilityCapacity[k]*yvariables["y_" +str(k)]
    prob +=s >= total_pop


    for i in ulist:
        s=""
        for j in centers:
            if i not in covered_units[j]: continue
            #if nodedij[i][j]>1.2*max_radius: continue
            s+=yvariables["y_" + str(j)]
        s -= zvariables["z_" +str(i)]
        prob += s>= 0

    if cover_pecentage >=99: 
        for i in ulist:
            prob += zvariables["z_" +str(i)] == 1
    else:
        s=""
        for i in ulist:
            noise=1.0
            if sampling==0:
                noise+=(random.random()-0.5)/100
            s+=nodes[i][3]*noise*zvariables["z_" +str(i)]
        prob += s>=total_pop*cover_pecentage/100.0
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message, timeLimit=time_limit,options=[("MIPGap",gap),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    if prob.status<=0:
        print "model unsolved..."
        return 0
    centersID=[-1 for x in range(num_districts)]
    if len(node_groups)<1: node_groups=[-1 for x in range(num_units)]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k] >=0:
                node_groups[i]=k
                break
    update_district_info()
    #for x in district_info: 
    #    if x[0]>0: print x
    return 1

def location_drop_cflpr(max_radius,cover_pecentage): 
    global centersID
    global node_groups
    clist=[x for x in range(num_districts) if centersID[x]>=0]
    if clist==[]: return 0
    dlist=[]
    random.shuffle(clist)
    for k in clist:
        ids=centersID[:]
        ids[k]=-1
        ulist=[[x,nodes[x][3]] for x in range(num_units) if node_groups[x]==k]
        ulist.sort(key=lambda x:-x[1])
        savings=district_info[k][2]
        savings+=facilityCost[k]
        caplist=[facilityCapacity[x]-district_info[x][1] for x in range(num_districts)]
        caplist[k]=0
        assigned=0
        for x,y in ulist:
            assigned=0
            for r in NearFacilityList[x]:
                if centersID[r]<0: continue
                if r==k: continue
                if y<=caplist[r]:
                    savings-=nodedik[x][r]
                    caplist[r]-=y
                    assigned=1
                    break
            if assigned==0: break
        if assigned==1: dlist.append([k,savings])
    if len(dlist)==0: return -1
    dlist.sort(key=lambda x:-x[1])
    update_cflpr_district_info(max_radius,cover_pecentage)
    f_del=0
    tabu=[]
    status=0
    while 1:
        if len(dlist)<1: break
        temp=[centersID[:],node_groups[:]]
        r=random.random()
        idx=int(r*0.999*min(3,len(dlist)))
        k=dlist[idx][0]
        del dlist[idx]
        if k in tabu: continue
        centersID[k]=-1
        ulist=[[x,nodes[x][3]] for x in range(num_units) if node_groups[x]==k]
        ulist.sort(key=lambda x:-x[1])
        caplist=[district_info[x][3]-district_info[x][1] for x in range(num_districts)]
        caplist[k]=0
        assigned=1
        for x,y in ulist:
            assigned=0
            for r in NearFacilityList[x]:
                if centersID[r]<0: continue
                if y<=caplist[r]:
                    caplist[r]-=y
                    node_groups[x]=r
                    assigned=1
                    if r not in tabu: tabu.append(r)
                    break
            if assigned==0: break
        if assigned==0:
            centersID=temp[0][:]
            node_groups=temp[1][:]
            update_cflpr_district_info(max_radius,cover_pecentage)
            return status
        obj=biobjective
        status=1
        update_cflpr_district_info(max_radius,cover_pecentage)
        f_del+=1
        if f_del>=len(clist)/20 or f_del>=10: break
        supply=sum(facilityCapacity[x] for x in range(num_districts) if centersID[x]>=0)
        if supply<=total_pop*3: break
    return status
def drop_method_cflpr(max_radius,cover_pecentage):
    global node_groups
    global centersID
    global all_solutions
    sol=[]
    node_groups=[NearFacilityList[x][0] for x in range(num_units)]
    centersID=facilityCandidate[:]
    for x in range(num_districts):
        if x not in node_groups: centersID[x]=-1
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k]==-1: continue
            node_groups[i]=k
            break 
    update_cflpr_district_info(max_radius,cover_pecentage)
    update_cflpr_centers(max_radius,cover_pecentage)
    update_cflpr_district_info(max_radius,cover_pecentage)
    while objective_overload>0:
        one_unit_move_cflpr(max_radius,cover_pecentage)
    while 1:
        sol=[centersID[:],node_groups[:]]
        sta=location_drop_cflpr(max_radius,cover_pecentage)  #too slow for large instances, sampling ?
        update_cflpr_district_info(max_radius,cover_pecentage)
        if sta<=0: break
        coverdemand=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
        cp=coverdemand*100.0/total_pop
        #if adaptive_number_of_facilities==1 and cp<cover_pecentage: break
        #print biobjective, objective_overload,cp
        if cp<cover_pecentage: 
            #centersID=sol[0][:]
            #node_groups=sol[1][:]
            #update_cflpr_district_info(max_radius,cover_pecentage)
            break
        sol=[biobjective,centersID[:],node_groups[:]]
        nf=sum(1 for x in centersID if x>=0)
        if adaptive_number_of_facilities==0 and nf <= max_num_facility: 
            sol=[biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0]
            break
    #if len(sol)>0:
    #    node_groups=sol[2][:]
    #    centersID=sol[1][:]
    #    update_cflpr_district_info(max_radius,cover_pecentage)
    nf=sum(1 for x in centersID if x>=0)
    if adaptive_number_of_facilities==0 and nf > max_num_facility:
        droplist=[]
        dlist=[x for x in range(num_districts) if centersID[x] >=0]
        while len(dlist) >max_num_facility:
            random.shuffle(dlist)
            k=dlist[0]
            dlist.remove(k)
            droplist.append(k)
            #print [len(dlist),len(droplist)],
        #print droplist
        for x in droplist: centersID[x]=-1
        ulist= [x for x in range(num_units) if node_groups[x] in droplist]
        for i in ulist:
            for k in NearFacilityList[i]:
                if centersID[k]<0: continue
                node_groups[i]=k
                break
        update_cflpr_district_info(max_radius,cover_pecentage)
        coverdemand=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
        cp=coverdemand*100.0/total_pop
        #print droplist,cp,biobjective,objective_overload
    update_cflpr_district_info(max_radius,cover_pecentage)
    return 1

def update_cflpr_district_info(max_radius,cover_pecentage):
    global biobjective
    update_district_info() 
    coverdemand=0
    for i in range(num_units):
        k=node_groups[i]
        if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
    if coverdemand  >= total_pop*cover_pecentage/100.0:
        return
    #obj2=(cover_pecentage-coverdemand*100.0/total_pop)*100000000*100
    obj1=pop_dis_coeff*objective_overload
    obj2=pop_dis_coeff*(total_pop*cover_pecentage/100.0-coverdemand)
    biobjective=biobjective+obj1+obj2

def update_cflpr_centers(r,p):
    global node_groups
    global centersID
    global time_update_centers
    t=time.time()
    centers=[x for x in range(num_districts) if centersID[x]>=0]
    random.shuffle(centers)
    for k in centers:
        newk=update_cflpr_center(k,r,p)
        if newk==k: continue
        ulist=[x for x in range(num_units) if node_groups[x]==k]
        for x in ulist: node_groups[x]=newk
        centersID[k]=-1
        centersID[newk]=facilityCandidate[newk]
        #print "a center updated!"
    time_update_centers+=time.time()-t

def update_cflpr_center(k,r,p):
    ulist=[x for x in range(num_units) if node_groups[x]==k]
    covered=0
    for i in ulist:
        j=node_groups[i]
        if nodedij[i][k]<r: 
            covered+=nodes[i][3]
    demand=sum(nodes[x][3] for x in ulist)
    best_cost=sum(nodedik[x][k] for x in ulist)
    best_center=k
    uid=nearCustomer[k]
    centers=[]
    for x in NearFacilityList[uid]:
        if centersID[x]<0: centers.append(x)
        if len(centers)>=min(20,len(ulist)): break
    for i in centers:
        if i==k: continue
        if facilityCapacity[i]<demand and facilityCapacity[k]<facilityCapacity[i]: continue
        cost=sum(nodedik[x][i] for x in ulist)
        if cost>=best_cost: continue
        covered3=0
        for j in ulist:
            if nodedij[j][i]<r:covered3+=nodes[j][3]
        if covered3<covered: continue 
        best_cost=cost
        best_center=i
    return best_center

#reduce overload only
def one_unit_move_cflpr(radius,c_percent): #FA
    global node_groups
    global time_op
    #if objective_overload>=0: return 0
    t=time.time()
    klist=[x for x in range(num_districts) if district_info[x][4]>0]
    ulist=[x for x in range(num_units) if node_groups[x] in klist]
    if ulist==[]:
        ulist=[x for x in range(num_units) if nodedij[x][node_groups[x]]>=radius]
    #if ulist==[]: return 0
    random.shuffle(ulist)
    ulist+=[x for x in range(num_units)]
    covered=0
    for x in range(num_units):
        k=node_groups[x]
        if nodedij[x][k]<=radius: covered+=nodes[x][3]
    max_coverd=total_pop*c_percent/100
    improved=0
    for nid in ulist:
        rid=node_groups[nid]
        demand=nodes[nid][3]
        klist=[x for x in NearFacilityList[nid] if centersID[x]>=0]
        klist=klist[:8]
        for k in klist:
            new_covered=covered
            if nodedij[nid][rid]<=radius: new_covered-=demand
            if nodedij[nid][k]<=radius: new_covered+=demand
            if new_covered<max_coverd and new_covered<covered: continue
            overload1=district_info[k][4]+district_info[rid][4]
            overload2=max(0, district_info[rid][1]-demand-facilityCapacity[rid])
            overload2+=max(0, district_info[k][1]+demand-facilityCapacity[k])
            if overload1 < overload2: continue
            if overload1 == overload2 and nodedij[nid][k]>=nodedij[nid][rid]: continue
            #if overload2==overload1: print nodedij[nid][k]-nodedij[nid][rid],
            node_groups[nid]=k
            update_cflpr_district_info(radius,c_percent)
            covered=new_covered
            #print int(objective),covered, objective_overload
            improved=1
            break
    update_cflpr_district_info(radius,c_percent)
    time_op[0]+=time.time()-t
    return improved

def cflpr_matheuristic(max_radius,cover_pecentage,multi_start,maxloops):
    global best_objective
    global best_biobjective
    global best_objective_overload
    global best_biobjective_global
    global best_centersID_global
    global best_solution_global
    global objective
    global biobjective
    global objective_overload
    global node_groups
    global centersID
    global district_info
    global facilityCost
    #global node_neighbors
    global region_pool
    global pool_index
    global is_spp_modeling
    global all_solutions
    global max_loops_solution_not_improved
    global multi_start_count
    global solver_message
    global pop_dis_coeff
    initialize_instance()
    max_loops_solution_not_improved=maxloops
    multi_start_count=multi_start
    print "num_facility:",max_num_facility
    pop_dis_coeff=100000000*1000.0/total_pop
    covered_demand_nodes=[]
    opt_radius=[0.0 for x in range(num_districts)]
    for k in range(num_districts):
        supply=facilityCapacity[k]
        dis=0.0
        covered_nodes=[]
        for i in nearCustomers[k]:
            if nodes[i][3]<supply:
                supply-=nodes[i][3]
                covered_nodes.append(i)
            else:
                dis=nodedij[i][k]
                covered_nodes.append(i)
                break
        covered_demand_nodes.append(covered_nodes[:])
        opt_radius[k]=dis
    #print "opt_radius[]",opt_radius
    varying_radius_max=[0.0 for x in range(num_units)]
    for i in range(num_units):
        maxd=0.0
        for k in range(num_districts):
            if i not in covered_demand_nodes[k]: continue
            if nodedij[i][k]>maxd: maxd=nodedij[i][k]
        #if maxd>max_radius: maxd=max_radius*5
        varying_radius_max[i]=maxd
    #ttt=[x for x in varying_radius_max if x<max_radius]
    #print varying_radius_max
    #print "varying_radius_max",len(ttt), sum(ttt)/len(ttt)
    nf=total_supply/total_pop/2
    for i in range(num_units):
        k=NearFacilityList[i][nf]
        if nodedij[i][k]>varying_radius_max[i]:
            varying_radius_max[i]=nodedij[i][k]
    #ttt=[x for x in varying_radius_max if x<max_radius]
    #print "varying_radius_max",len(ttt), sum(ttt)/len(ttt)

    #print "varying_radius_max",
    #for i in range(num_units):
    #    print i,varying_radius_max[i]
    all_solutions=[]
    region_pool=[]
    t=time.time()
    best_biobjective_global = MAXNUMBER
    best_biobjective = MAXNUMBER
    district_info = [[0,0.0,0,0,0] for x in range(num_districts)]
    population=[] #all
    pool_index=[[] for x in range(num_districts)]
    for x in range(num_districts):
        if facilityCost[x]<1000000:  facilityCost[x]=100000000
    not_improved_global=0
    init_methods=[0,1,2] #SCLP, CSCLP, Drop
    if adaptive_number_of_facilities==0: init_methods=[1,2]
    loops=multi_start_count
    if initial_solution_method<0 and multi_start_count<6:
        loops=6
    for idx in range(loops):
        #cflpr_lp_model(max_radius,cover_pecentage,-1,100,0.000000000001)
        #add_method_cflpr2(max_radius,cover_pecentage)
        #CFLP_CKFLP_model_init(max_radius,500,0.01, 1)
        #lscp_mip_model_init(max_radius,cover_pecentage,100,0.01)
        #drop_method_cflpr(max_radius,cover_pecentage)
        init_t=time.time()
        if initial_solution_method<0:
            nm=len(init_methods)
            method=init_methods[idx%nm]
        else:
            method=initial_solution_method
        if method==0: 
            sta=lscp_mip_model_init(max_radius,cover_pecentage,100,0.002)
            if sta<=0: 
                if initial_solution_method<0: init_methods.remove(0)
                continue
        if method==1: 
            sta=clscp_mip_model_init(max_radius,cover_pecentage,100,0.002)
            if sta<=0: 
                if initial_solution_method<0: init_methods.remove(1)
                continue
        if method==2: 
            sta=drop_method_cflpr(max_radius,cover_pecentage)
        update_cflpr_district_info(max_radius,cover_pecentage)
        #print biobjective,
        #update_cflpr_centers(max_radius,cover_pecentage)
        #update_cflpr_district_info(max_radius,cover_pecentage)
        #print biobjective,
        #print "time",int(time.time()-init_t),
        while one_unit_move_cflpr(max_radius,cover_pecentage):
            pass
        #print biobjective,
        #print int(time.time()-init_t),
        coverdemand=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
        cp=coverdemand*100.0/total_pop
        nf=sum(1 for x in centersID if x>=0)
        update_best_solution()
        if initial_solution_method<0 and adaptive_number_of_facilities==1 and objective_overload>total_pop/20:
            if method==0 and method in init_methods and len(init_methods)>1: 
                 init_methods.remove(method)
            elif biobjective>best_biobjective_global*1.5 and len(init_methods)>1: 
                 init_methods.remove(method)
        print "init. sol.",idx, "mthd",method,nf,biobjective,objective,objective_overload, int(cp*100)/100.0,time.time()-t
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0])
    population.sort(key=lambda x:x[0]+x[-1]*10000)
    best_biobjective_global=population[0][0]
    best_centersID_global=population[0][1][:]
    best_solution_global=population[0][2][:]
    #best_biobjective_global = MAXNUMBER
    #best_biobjective = MAXNUMBER
    update_cflpr_district_info(max_radius,cover_pecentage)
    update_best_solution()
    all_solutions=population
    loop=-1
    while 1:
        r=random.random()
        sidx = int(min(multi_start_count-1,len(population))* pow(r,1.5)*0.999)
        node_groups=population[sidx][2][:]
        centersID=population[sidx][1][:]
        update_cflpr_district_info(max_radius,cover_pecentage)
        loop+=1
        not_improved_global+=1
        obj=best_biobjective_global
        t_mip=time.time()
        r=random.random()
        if r<0.33:
            sscflpr_sub_assign_model(max_radius,cover_pecentage,20,0.0000000001)
        else:
            #cflpr_sub_model(max_radius,cover_pecentage,20,0.0000000001)
            sscflpr_sub_model(max_radius,cover_pecentage,20,0.0000000001) #20s?? or 50s?
        update_cflpr_district_info(max_radius,cover_pecentage)
        update_best_solution()
        update_region_pool_all()
        print int(time.time()-t_mip),
        coverdemand=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
        cover=min(0,coverdemand-total_pop*cover_pecentage/100)
        overload=objective_overload
        while one_unit_move_cflpr(max_radius,cover_pecentage):
            pass
        coverdemand=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
        print "M",min(0,coverdemand-total_pop*cover_pecentage/100)-cover,int(objective_overload-overload),
        update_best_solution()
        update_region_pool_all()

        if objective_overload>=10*10000000000:
            temp=[best_biobjective_global,best_centersID_global[:],best_solution_global[:]]
            VND_local_search()
            update_cflpr_district_info(max_radius,cover_pecentage)
            update_best_solution()
            update_region_pool_all()
            if biobjective >temp[0]:
                best_biobjective_global=temp[0]
                best_centersID_global=temp[1][:]
                best_solution_global=temp[2][:]


        temp=[best_biobjective_global,best_centersID_global[:],best_solution_global[:]]
        update_cflpr_centers(max_radius,cover_pecentage)
        update_cflpr_district_info(max_radius,cover_pecentage)
        update_best_solution()
        update_region_pool_all()
        if biobjective >temp[0]:
            best_biobjective_global=temp[0]
            best_centersID_global=temp[1][:]
            best_solution_global=temp[2][:]

        #print int(objective),
        #update_cflpr_centers()
        #update_cflpr_district_info(max_radius,cover_pecentage)
        #update_best_solution()

#        temp=[best_biobjective_global,best_centersID_global[:],best_solution_global[:]]
#        assign_ruin_recreate(-1)
#        VND_local_search()
#        update_cflpr_district_info(max_radius,cover_pecentage)
#        update_best_solution()
#        update_region_pool_all()
#        if adaptive_number_of_facilities==0:
#            nf=sum(1 for x in centersID if x>=0)
#            if nf<max_num_facility: continue

        #population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,max(0,total_pop-objective_supply)])
        s="" 
        if best_biobjective_global<obj-0.000001:  #0.1%
            s="*"
            impp=((best_biobjective_global-obj)/best_biobjective_global)*1000 #-0.1%
            not_improved_global+=int(max_loops_solution_not_improved*impp)
            if not_improved_global<0: not_improved_global=0
        else: s="-"
        bnf=sum(1 for x in best_centersID_global if x>=0)
        s+="Loop "+str(loop) + " Best: " +str(bnf)+" "+str(int(best_biobjective_global))+" "+str(int(best_objective_global))+" "+str(int(best_overload_global))
        bnf=sum(1 for x in centersID if x>=0)
        s+=" Current: "+str(bnf)+" "+str(int(biobjective))+" "+str(int(objective))+" "+str(int(objective_overload))
        s+=" Info: " +str(not_improved_global)+ " "+str(int(time.time()-t))
        coverdemand=0
        for i in range(num_units):
            k=node_groups[i]
            if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
        cp=coverdemand*1.0/total_pop
        s+=" "+str(int(cp*10000)/100.0)
        arcpy_print(s)
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0])
        population.sort(key=lambda x:x[0])#+x[-1]*10000)
        population=pop_selection(population)
        all_solutions=population
        if not_improved_global >= max_loops_solution_not_improved: break
        #if time.time()-t_ga > heuristic_time_limit:  break
    #post procedure
    ga_time=time.time()-t
    print "Heuristic solution:",biobjective,objective_fcost,objective,objective_overload,ga_time
    t_spp=time.time()
    if is_spp_modeling>=1:
        arcpy_print("SPP modelling..."+str(len(region_pool)) )
        lspp_model(max_radius,cover_pecentage,ga_time/10,0.00000001)
        update_cflpr_district_info(max_radius,cover_pecentage)
        while one_unit_move_cflpr(max_radius,cover_pecentage):
            pass
        update_best_solution()
        print "spp solution:",biobjective,objective_fcost,objective,objective_overload,time.time()-t_spp
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,max(0,total_pop-objective_supply)])
    population.sort(key=lambda x:x[0]+x[5]*1000000)
    #population.sort(key=lambda x:x[0])
    coverdemand=0
    for i in range(num_units):
        k=node_groups[i]
        if nodedij[i][k]<=max_radius: coverdemand+=nodes[i][3]
    node_groups=best_solution_global[:]
    centersID=best_centersID_global[:]

    print "final solution:",biobjective,objective_fcost,objective,objective_overload,objective/total_pop,int(coverdemand*10000/total_pop)/100.0
    return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]

#k_medoids() #max_num_facility
#
def sscflpr_sub_model(max_radius,cover_pecentage,time_limit,mipgap): 
    #may be more facility needed, due to the covering constraint 
    global centersID
    global node_groups
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    u=-1
    if objective_overload>0:
        rlist=[x for x in range(num_districts) if centersID[x]>=0 and district_info[x][4]>0]
        ulist=[x for x in range(num_units) if node_groups[x] in rlist]
        random.shuffle(ulist)
        u=ulist[0]
    dlist=select_region(-1) #select_region_cflpr(-1)#select_region(-1) 
    ulist=[x for x in range(num_units) if node_groups[x] in dlist]
    cloclist=get_cand_locations(dlist,ulist)
    overl=sum(district_info[x][4] for x in dlist)
    #print overl,
    centers=[]
    for clist in cloclist:
        tblist=list(set(clist+dlist))
        tblist.sort()
        centers=tblist
        break
    #print [len(centers),len(ulist)],
    uid=-1
    udis=MAXNUMBER
    for x in ulist:
        dis=sum(nodedij[x][k] for k in dlist)
        if dis < udis:
            udis=dis
            uid=x
    exist_facility=[] #x for x in range(num_districts) if centersID[x]>=0 and x not in dlist and facilityCapacity[x] > district_info[x][1]]

    for k in NearFacilityList[uid]:
        if centersID[k]<0: continue
        if district_info[k][1]>=facilityCapacity[k]: continue
        if k in dlist: continue
        exist_facility.append(k)
        if len(exist_facility)>=len(dlist): break
    prob = pulp.LpProblem("pcp",pulp.LpMinimize)
    xvariables={}
    costs={}
    for i in ulist:
        for j in centers:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]#*nodedij[i][j]
        for j in exist_facility:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]

    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        costs["y_" +str(i)]=100000000 #facilityCost[i]
    zvariable=pulp.LpVariable("z", 0, None, pulp.LpContinuous)
    overload=max(0,objective_overload)
    z2variables={}
    for i in centers:
        z2variables["z2_" +str(i)]=pulp.LpVariable("z2_" +str(i), 0, None, pulp.LpInteger)
    obj=""
    for x in yvariables:
        obj +=costs[x]* yvariables[x]
    for i in ulist:
         for j in centers:
            obj+=costs["x_" +str(i)+ "_"+ str(j)]*xvariables["x_" +str(i)+ "_"+ str(j)]
         for j in exist_facility:
            obj+=costs["x_" +str(i)+ "_"+ str(j)]*xvariables["x_" +str(i)+ "_"+ str(j)]
    penalty=pop_dis_coeff #100000000*10000.0/total_pop #pop_dis_coeff #pop_dis_coeff*objective/total_pop
    if cover_pecentage>=99: penalty*=10
    for x in z2variables:
        obj += penalty*z2variables[x]
    obj+=penalty*zvariable #100*100000000*zvariable #*sum(facilityCost)/num_districts*zvariable
    prob += obj

    for i in ulist:
        s=""
        for j in centers:
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        for j in exist_facility:
            v="x_" +str(i)+ "_"+ str(j)
            if v in xvariables : s+=xvariables[v]
        prob += s==1

    for k in centers:
        s=""
        for i in ulist:
            s+= nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(k) ]
        s-= facilityCapacity[k]*yvariables["y_" +str(k)]
        s-= z2variables["z2_" +str(k)]
        prob += s <= 0

    for k in centers:
        for i in ulist:
            s= z2variables["z2_" +str(k)]
            s-= facilityCapacity[k]*yvariables["y_" +str(k)]
            prob += s <= 0

    for k in exist_facility:
        s=""
        for i in ulist:
            v="x_" +str(i)+ "_"+ str(k)
            if v in xvariables: s+=nodes[i][3]*xvariables[v]
        prob+= s <= facilityCapacity[k] - district_info[k][1]

    total_cover1=0
    total_cover2=0
    for i in range(num_units):
        if nodedij[i][node_groups[i]]>max_radius: continue
        if i in ulist: 
           total_cover1+=nodes[i][3]
        else: total_cover2+=nodes[i][3]
    min_covered=total_pop*cover_pecentage/100
    if total_cover2<min_covered:
        s=""
        for i in ulist:
            for j in centers:
                if nodedij[i][j]>max_radius: continue
                s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
            for j in exist_facility:
                if nodedij[i][j]>max_radius: continue
                s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
        s+=zvariable
        prob += s>=min_covered-total_cover2

    s=""
    for x in yvariables: s+=yvariables[x]
    numf=sum(1 for x in centersID if x>=0)
    if adaptive_number_of_facilities==0:
        if numf==max_num_facility:
            prob += s== len(dlist)
    else:
        ngap=max(1,nf*5/100)
        prob += s>= len(dlist)-ngap
        prob += s<= len(dlist)+ngap
    #prob.writeLP("_SSCFLP.lp")
    initvalues=loc_sub_mst(dlist,ulist)
    for x,v in initvalues:
        if x.find("x")==0:  xvariables[x].setInitialValue(v)
        if x.find("y")==0:  yvariables[x].setInitialValue(v)
    #warmStart=True,
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap),'set mip tolerances absmipgap 100', 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit,options=[("MIPGap",gap),("MIPGapAbs",gap*objective),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    if prob.status<=0:
        print "model unsolved... or infeasible, msg=", prob.status
        prob.writeLP("_sscflpr_sub_model.lp")
        return []
    nf=0
    nlist=[]
    obj1=0.0
    for x in ulist: obj1+=nodedik[x][node_groups[x]]
    for x in dlist: centersID[x]=-1
    for x in ulist: node_groups[x]=-1
    for v in prob.variables():
        if (v.varValue >= 0.9):
            if v.name=="z": continue
            if v.name=="z2": continue
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
                nlist.append(i)
                nf+=1
            if items[0]=='x':
                node_groups[i]=int(items[2])
    #update_district_info()
    update_cflpr_district_info(max_radius,cover_pecentage)
    obj2=0.0
    for x in ulist: obj2+=nodedik[x][node_groups[x]]
    print len(dlist),len(centers),nf, int(obj1),"->",int(obj2),
    #if nf<-100: 
    #    dlist.sort()
    #    nlist.sort()
    #    centers.sort()
    #    print pulp.value(prob.objective),dlist,nlist,centers
    #    print initvalues
    return 1
def loc_sub_mst(dlist,ulist):
    vars=[]
    for i in ulist:
        k=node_groups[i]
        if k<0: continue
        v='x_'+str(i)+ '_'+ str(k)
        vars.append([v,1])
    for i in dlist:
        if centersID[i]<0:continue
        v='y_'+str(i)
        vars.append([v,1])
    return vars

def sscflpr_sub_assign_model(max_radius,cover_pecentage,time_limit,mipgap): 
    #may be more facility needed, due to the covering constraint 
    global centersID
    global node_groups
    uid=random.randint(0,num_units-1)
    centers=[] #x for x in range(num_districts) if centersID[x]>=0 and x not in dlist and facilityCapacity[x] > district_info[x][1]]
    for k in NearFacilityList[uid]:
        if centersID[k]<0: continue
        centers.append(k)
        if len(centers)>=min(16,num_districts/10):
             break
    ulist=[x for x in range(num_units) if node_groups[x] in centers]

    prob = pulp.LpProblem("gap",pulp.LpMinimize)
    xvariables={}
    costs={}
    for i in ulist:
        for j in centers:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]#*nodedij[i][j]
    z=pulp.LpVariable("z", 0, None, pulp.LpInteger)
    zvariables={}
    for j in centers:
        zvariables["z_" + str(j)]=pulp.LpVariable("z_" + str(j), 0, None, pulp.LpInteger)
    obj=""
    for i in ulist:
         for j in centers:
            obj+=costs["x_" +str(i)+ "_"+ str(j)]*xvariables["x_" +str(i)+ "_"+ str(j)]
    penalty=pop_dis_coeff*objective/total_pop
    for j in centers:
        obj+=penalty*zvariables["z_" + str(j)]
    obj+=penalty*z
    prob += obj

    for i in ulist:
        s=""
        for j in centers:
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        prob += s==1

    for k in centers:
        s=""
        for i in ulist:
            s+= nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(k) ]
        s-=zvariables["z_" + str(k)]
        prob += s <= facilityCapacity[k]
    total_cover1=0
    total_cover2=0
    for i in range(num_units):
        if nodedij[i][node_groups[i]]>max_radius: continue
        if i in ulist: 
           total_cover1+=nodes[i][3]
        else: total_cover2+=nodes[i][3]
    s=""
    for i in ulist:
        for j in centers:
            if nodedij[i][j]>max_radius: continue
            s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
    min_covered=total_pop*cover_pecentage/100
    s+=z
    prob += s>=min_covered-total_cover2

    #max_covered=total_pop*(cover_pecentage-0.1)/100
    #if total_cover1+total_cover2 >= max_covered:
    #    prob += s>=max_covered-total_cover2
    #if total_cover1+total_cover2 < max_covered:
    #    prob += s>=total_cover1

    initvalues=loc_sub_mst(centers,ulist)
    for x,v in initvalues:
        if x.find("x")==0:  xvariables[x].setInitialValue(v)
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap),'set mip tolerances absmipgap 100', 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit,options=[("MIPGap",gap),("MIPGapAbs",gap*objective),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    if prob.status<=0:
        print "model unsolved... or infeasible, msg=", prob.status
        prob.writeLP("_sscflpr_sub_assign.lp")
        return []
    nf=len(centers)
    nlist=[]
    obj1=0.0
    for x in ulist: obj1+=nodedik[x][node_groups[x]]
    for x in ulist: node_groups[x]=-1
    for v in prob.variables():
        if (v.varValue >= 0.9):
            items=v.name.split('_')
            if items[0]=='x':
                i=int(items[1])
                node_groups[i]=int(items[2])
    obj2=0.0
    for x in ulist: obj2+=nodedik[x][node_groups[x]]
    print len(centers),len(centers),nf, int(obj1),"->",int(obj2),
    return 1

def cflpr_sub_model(max_radius,cover_pecentage,time_limit,mipgap): 
    #may be more facility needed, due to the covering constraint 
    global centersID
    global node_groups
    dlist=select_region_cflpr(-1)#select_region(-1) 
    ulist=[x for x in range(num_units) if node_groups[x] in dlist]
    cloclist=get_cand_locations(dlist,ulist)
    centers=[]
    for clist in cloclist:
        tblist=list(set(clist+dlist))
        tblist.sort()
        centers=tblist
        break
    #print [len(centers),len(ulist)],
    uid=-1
    udis=MAXNUMBER
    for x in ulist:
        dis=sum(nodedij[x][k] for k in dlist)
        if dis < udis:
            udis=dis
            uid=x
    exist_facility=[] #x for x in range(num_districts) if centersID[x]>=0 and x not in dlist and facilityCapacity[x] > district_info[x][1]]
    for k in NearFacilityList[uid]:
        if centersID[k]<0: continue
        if k in dlist: continue
        exist_facility.append(k)
        if len(exist_facility)>=len(dlist): break
    prob = pulp.LpProblem("pcp",pulp.LpMinimize)
    xvariables={}
    costs={}
    for i in ulist:
        for j in centers:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpContinuous)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]#*nodedij[i][j]
        for j in exist_facility:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpContinuous)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]
    yvariables={}
    for i in centers:
        yvariables["y_" +str(i)]=pulp.LpVariable("y_" +str(i), 0, 1, pulp.LpBinary)
        costs["y_" +str(i)]=100000000#facilityCost[i]
    zvariable=pulp.LpVariable("z", 0, None, pulp.LpContinuous)

    obj=""
    for x in yvariables:
        obj +=costs[x]* yvariables[x]
    #obj+=100000*dvariable
    for i in ulist:
         for j in centers:
            obj+=costs["x_" +str(i)+ "_"+ str(j)]*xvariables["x_" +str(i)+ "_"+ str(j)]
         for j in exist_facility:
            obj+=costs["x_" +str(i)+ "_"+ str(j)]*xvariables["x_" +str(i)+ "_"+ str(j)]
    penalty=pop_dis_coeff#*objective/total_pop
    obj+=penalty*zvariable
    prob += obj

    for i in ulist:
        s=""
        for j in centers:
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        for j in exist_facility:
            v="x_" +str(i)+ "_"+ str(j)
            if v in xvariables : s+=xvariables[v]
        prob += s==1

    for k in centers:
        s=""
        for i in ulist:
            s+= nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(k) ]
        s-= facilityCapacity[k]*yvariables["y_" +str(k)]
        prob += s <= 0
    for k in exist_facility:
        s=""
        for i in ulist:
            v="x_" +str(i)+ "_"+ str(k)
            if v in xvariables: s+=nodes[i][3]*xvariables[v]
        prob+= s <= facilityCapacity[k] - district_info[k][1]
    total_cover1=0
    total_cover2=0
    for i in range(num_units):
        if nodedij[i][node_groups[i]]>max_radius: continue
        if i in ulist: 
           total_cover1+=nodes[i][3]
        else: total_cover2+=nodes[i][3]
    s=""
    for i in ulist:
        for j in centers:
            if nodedij[i][j]>max_radius: continue
            s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
        for j in exist_facility:
            if nodedij[i][j]>max_radius: continue
            s+=nodes[i][3]*xvariables["x_" +str(i) + "_"+ str(j)]
    s+=zvariable
    prob += s >= total_pop*cover_pecentage/100-total_cover2
    #prob += s>=min(total_cover1,max(0,total_pop*cover_pecentage/100-total_cover2))
    #if total_cover1+total_cover2 >=total_pop*cover_pecentage/100: 

    if adaptive_number_of_facilities==0:
        numf=sum(1 for x in centersID if x>=0)
        if numf==max_num_facility:
            s=""
            for x in yvariables: s+=yvariables[x]
            prob += s== len(dlist)
    #prob.writeLP("_pcp.lp")
    initvalues=loc_sub_mst(dlist,ulist)
    for x,v in initvalues:
        if x.find("x")==0:  xvariables[x].setInitialValue(v)
        if x.find("y")==0:  yvariables[x].setInitialValue(v)
    #warmStart=True,
    gap=mipgap
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit, options=['set mip tolerances mipgap '+ str(gap),'set mip tolerances absmipgap 100', 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message,warmStart=True, timeLimit=time_limit,options=[("MIPGap",gap),("MIPGapAbs",gap*objective),("TimeLimit",time_limit)])
    solver.setTmpDir()
    solver.actualSolve(prob)
    if prob.status<=0:
        print "model unsolved... or infeasible, msg=", prob.status
        return []
    nf=len(dlist)
    nlist=[]
    xij=[0.0 for x in range(num_units)]
    obj1=0.0
    for x in ulist: obj1+=nodedik[x][node_groups[x]]
    for x in dlist: centersID[x]=-1
    for x in ulist: node_groups[x]=-1
    for v in prob.variables():
        if (v.varValue >= 0.000001):
            if v.name=="z": continue
            items=v.name.split('_')
            i=int(items[1])
            if items[0]=='y':
                centersID[i]=facilityCandidate[i]
                nlist.append(i)
                nf-=1
            if items[0]=='x':
                if v.varValue >xij[i]:
                    xij[i]=v.varValue
                    node_groups[i]=int(items[2])
    #update_district_info()
    update_cflpr_district_info(max_radius,cover_pecentage)
    obj2=0.0
    for x in ulist: obj2+=nodedik[x][node_groups[x]]
    print len(dlist),len(centers),nf, int(obj1),"->",int(obj2),
    #if nf<-100: 
    #    dlist.sort()
    #    nlist.sort()
    #    centers.sort()
    #    print pulp.value(prob.objective),dlist,nlist,centers
    #    print initvalues
    return 1

def lspp_model(maxr,cover_pecentage,maxtime,mipgap): #bug, the covering constraint
    global node_groups
    global centersID
    if len(region_pool)<=10:
        ##noprint "no candidate district!"
        print "len(region_pool)<=10", len(region_pool)
        return 0
    alpha_coeff=avg_dis_min*pop_dis_coeff  
    prob = pulp.LpProblem("sdp_spp",pulp.LpMinimize)
    variables=[]
    costs=[]
    for i in range(len(region_pool)):
        x=pulp.LpVariable("x_" +"{0:07d}".format(i), 0, 1,pulp.LpBinary)
        variables.append(x)
        cost=region_pool[i][1]+region_pool[i][3]*100000000#alpha_coeff
        k=region_pool[i][4]
        c2=sum(nodedik[x][k] for x in region_pool[i][0])
        if abs(c2-cost)>0.001: print "check...", i,k,c2,cost,region_pool[i][1],region_pool[i][3]
        if fixed_cost_obj==1: ##[ulist,dis,demand,overload,k]
            cost+=facilityCost[k]
        costs.append(cost)

    obj=""
    for i in range(len(region_pool)):
        obj+=costs[i]*variables[i]
    prob+=obj

    maxr_cover_list=[0 for i in range(len(region_pool))]
    rlist=[[] for i in range(num_units)]
    for j in range(len(region_pool)):
        k=region_pool[j][4]
        for x in region_pool[j][0]:
            rlist[x].append(j)
            if nodedij[x][k]>maxr: maxr_cover_list[j]+= nodes[x][3]
    for i in range(num_units):
        s=""
        for x in rlist[i]:
            s+=variables[x]
        if spatial_contiguity==0:
            prob+=s >= 1
        else:
            prob+=s == 1
    if spatial_contiguity==0:
        for k in range(num_districts):
            s=""
            for i in range(len(region_pool)):
                if region_pool[i][4]!=k: continue
                s+=variables[i]
            if len(s)>0: prob+=s <= 1
    s=""
    for i in range(len(region_pool)):
        if maxr_cover_list[i]>0:
            s+=maxr_cover_list[i]*variables[i]
    if len(s)>1:
        maxc=total_pop - total_pop*cover_pecentage/100
        prob += s<=maxc
        #print s, maxc
    if adaptive_number_of_facilities==0:
        s=""
        for i in range(len(variables)):
            s+=variables[i]
        prob+= s==max_num_facility
    #prob.writeLP("_spp.lp")
    #mip_mst_file=tempfile.mkstemp()[1].split("\\")[-1]

    vars=spp_mst()
    for x,v in vars: variables[x].setInitialValue(v)
    solver=0
    if mip_solver=='cbc': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.COIN_CMD(timeLimit=maxtime,mip=1,msg=solver_message,gapRel=mipgap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.cplex_api.CPLEX_CMD(msg=solver_message,warmStart=True,timeLimit=maxtime,options=['set parallel -1','set mip tolerances mipgap ' + str(mipgap)])
    if mip_solver=='gurobi': #solver_message #'set emphasis mip 3','set threads 4', 
        solver=pulp.apis.GUROBI_CMD(msg=solver_message,warmStart=True,options=[("MIPGap",mipgap),("TimeLimit",maxtime)])
    solver.setTmpDir() #=mip_file_path
    solver.actualSolve(prob)
    #if os.path.isfile(mip_mst_file): os.remove(mip_mst_file)
    if prob.status<=0:
       print "no solution! prob.status<0..."
       return prob.status #failer
    node_groups=[-1 for x in range(num_units)]
    centersID=[-1 for x in range(num_districts)]
    for v in prob.variables():
        if (v.varValue >= 0.99):
            items=v.name.split('_')
            i=int(items[1])
            k=region_pool[i][4]
            #print k,costs[i],facilityCost[k]
            centersID[k]=facilityCandidate[k]
            for x in region_pool[i][0]:
                node_groups[x]=k
    #a=[ x for x in centersID if x>=0]
    #print "spp locs:",len(a),a
    update_district_info()
    #for i in range(num_districts): 
    #    if district_info[i][1] >0: print i,district_info[i],facilityCost[i]
    #print 
    return 1 #success
def spp_mst():
    vars=[]        
    num_pool=len(region_pool)
    for k in range(num_districts):
        ulist=[x for x in range(num_units) if node_groups[x]==k] 
        for i in range(num_pool):
            if ulist==region_pool[i][0] and region_pool[i][4]==k:
                vars.append([i,1])
                break
    return vars

def pop_selection(population):
    population1=copy.deepcopy(population)
    population1.sort(key=lambda x:x[0])
    population2=[] #delete identical solution
    #sol=[best_biobjective_global,best_centersID_global[:],best_solution_global[:],best_objective_global,best_objective_fcost_global,best_overload_global,0]
    #population2.append(copy.deepcopy(sol))
    population2.append(copy.deepcopy(population1[0]))
    sdiff=1
    if location_problem==3:
        sdiff=max_num_facility*5/100
        if sdiff<=5: sdiff=5
        
    for x in population1:
        issimilar=0
        for y in population2:
            rlist=[i for i in range(num_districts) if x[1][i] != y[1][i]]
            if len(rlist)>=sdiff: continue
            else:
                if location_problem>=1:
                    issimilar=1
                    break
            ulist=[i for i in range(num_units) if x[2][i] != y[2][i]]
            #diffpop=sum(nodes[i][3] for i in ulist)
            #if len(ulist)<min(num_units*1.0/num_districts,num_units/30.0) and diffpop*100.0/total_pop < min(3.0,100.0/num_districts): #100.0/num_districts: #<5%
            #print len(ulist),
            if len(ulist)<num_units*(solution_similarity_limit/100.0):
                issimilar=1
                break
        if issimilar==0:
            population2.append(copy.deepcopy(x))
        if len(population2)>=max(multi_start_count*2,10):
            break
    return population2


def update_centers():
    global node_groups
    global centersID
    global time_update_centers
    if location_problem==1 or location_problem==0: return
    t=time.time()
    obj=biobjective
    sol=[-1 for x in range(num_units)]
    centers=[]
    for k in range(num_districts):
        if centersID[k]==-1: continue
        kn,ulist=update_center(k)
        for x in ulist: sol[x]=kn
        centers.append(kn)
        centersID[k]=-1
        #print [k,kn,k in ulist],
    node_groups=sol[:]
    for k in centers:
        centersID[k]=facilityCapacity[k]
    #if location_problem==3:
    #    for i in range(num_units):
    #        for k in NearFacilityList[i]:
    #            if centersID[k]>=0:
    #                node_groups[i]=k
    #                break
    obj=biobjective
    update_district_info()
    update_best_solution()
    #print "updatecenters",biobjective-obj
    time_update_centers+=time.time()-t
    #print obj,obj-biobjective

def update_center(k):
    if location_problem==3: #pmp
        return update_center_pmp(k)
    ulist=[x for x in range(num_units) if node_groups[x]==k]
    if ulist==[]: return k,[]
    best_cost=sum(nodedik[x][k] for x in ulist)
    best_center=k
    if pmp_I_eaquls_J==1:
        for i in ulist:
            cost=sum(nodedik[x][i] for x in ulist)
            if cost<best_cost:
                best_cost=cost
                best_center=i
    if pmp_I_eaquls_J==0:
        for i in range(num_districts):
            cost=sum(nodedik[x][i] for x in ulist)
            if cost<best_cost:
                best_cost=cost
                best_center=i
    return best_center,ulist


# return k centers by k-means
def k_medoids_cflpr_sampling(num_k):
    global centersID
    global node_groups
    centers=[]
    sol=[-1 for x in range(num_units)]
    distance_obj=0.0
    while 1:
        if len(centers)==num_k: break
        k=random.randint(0,num_districts-1)
        if k not in centers: centers.append(k)
    #for i in range(num_units):
    #    for k in NearFacilityList[i]:
    #        if k in centers:
    #            distance_obj+=nodedik[i][k]
    #            sol[i]=k
    #            break
    centers2=[]
    sol2=[-1 for x in range(num_units)]
    distance_obj2=0.0
    loop=0
    #k-means
    while 1:
        loop+=1
        distance_obj=0.0
        for i in range(num_units):
            for k in NearFacilityList[i]:
                if k in centers:
                    sol[i]=k
                    distance_obj+=nodedik[i][k]
                    break
        distance_obj2=0.0
        centers2=[]
        for k in centers:
            ulist=[x for x in range(num_units) if sol[x]==k]
            if len(ulist)==0: continue
            clist=range(num_districts)
            cid=-1
            mindis=MAXNUMBER
            for i in clist:
                dsum=sum(nodedik[x][i] for x in ulist)
                if dsum<mindis:
                    mindis=dsum
                    cid=i
            centers2.append(cid)
            distance_obj2+=mindis
            for i in ulist: sol2[i]=cid
        sol=sol2[:]
        centers=list(set(centers2))
        #print distance_obj,distance_obj2
        if abs(distance_obj2-distance_obj)/distance_obj<0.0001:
            break
    return centers

def location_check(key):
    if -1 in node_groups:
        arcpy_print("debug: "+str(key)+ " unit(s) not assigned! ")
        #return -1
    rlist=list(set(node_groups))
    if -1 in rlist: 
        arcpy_print("debug: "+str(key)+ " demand not assigned"+str(rlist))
        arcpy_print(str(node_groups))
        rlist.remove(-1)
    if len(rlist)>max_num_facility and adaptive_number_of_facilities==0:
        arcpy_print("debug: "+str(key)+ " too many facilities"+str(rlist))
        #return -1
    for k in range(num_districts):
        if k in rlist and centersID[k]==-1:
            arcpy_print("debug: "+str(key)+ " facilitiy not selected but used")
            #return -1
        if centersID[k]>=0 and k not in rlist:
            arcpy_print("debug: "+str(key)+ " facilitiy selected but not used")
            print k, district_info[k]
            print [x for x in centersID if x>=0]

            #return -1
        uid=centersID[k]
        if spatial_contiguity==1 and uid>-1 and node_groups[uid]!=k:
            arcpy_print("debug: "+str(key)+ " facilitiy unit assigned to other facility"+str(k))
            print k,uid, node_groups[uid]
            #return -1
    #return 1

def print_solution():
    arcpy_print("_______________final solution_______________")
    for i in range(num_units):
        s=""
        for x in nodes[i]:
            s+=str(x)+"\t"
        k=node_groups[i]
        kunit=centersID[k]
        s+=str(nodes[kunit][4])+"\t"
        selected=-1
        if i in facilityCandidate:
            selected=0
            if i in centersID:
                selected=1
        s+=str(selected)
        arcpy_print(s)

def search_stat():
    arcpy_print("----------------search statistics----------------------")
    arcpy_print("one unit move, move and time: "+ str(count_op[0])+ ", " +str(time_op[0]) )
    arcpy_print("two unit move, move and time: "+ str(count_op[1])+ ", " +str(time_op[1]) )
    arcpy_print("three unit move, move and time: "+ str(count_op[2])+ ", " +str(time_op[2]) )
    arcpy_print("location swap time: "+ str(time_location[0]) )
    arcpy_print("location drop time: "+ str(time_location[1]) )
    arcpy_print("location add time: "+ str(time_location[2]) )
    arcpy_print("location add-drop time: "+ str(time_location[3]) )
    arcpy_print("location multi-exchange time: "+ str(time_location[4]) )
    arcpy_print("r_r_reselect_location_pmp time: "+ str(time_location[5]) )
    arcpy_print("location TB heur. time: "+ str(time_location[7]) )
    arcpy_print("TB_Whitaker time: "+str(time_Whitaker))
    arcpy_print("location PMP sub_mip time: "+ str(time_location[8]) )
    arcpy_print("location CFLP sub_mip time: "+ str(time_location[9]) )
    arcpy_print("location PMP TB time: "+ str(time_pmp_re_location) )
    arcpy_print("repair time: "+ str(time_repair) )
    arcpy_print("check edge unit time: "+str(time_check_edge_unit))
    arcpy_print("update_centers time: "+ str(time_update_centers) )
    arcpy_print("spp regions: "+ str(len(region_pool)) )
    arcpy_print("spp pooling time: "+ str(time_spp) )
    arcpy_print("connectivity check time: "+ str(time_check))
    arcpy_print("time for ruin and recreate: " + str(time_ruin_recreate))
    if spatial_contiguity==1:
        sta=check_solution_continuality_feasibility(best_solution_global)
        arcpy_print("solution on continuality (0 no, 1 yes) : "+str(sta))
    arcpy_print("----------------end of search statistics----------------")
