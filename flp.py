import SSCFLPr as d
import time,sys,random
maxr=1.0
maxcp=80
spp=0 #set partitioning
psize=10 #population size
maxloops=100 #max_loops_solution_not_improved
t=50 #timelimit for searching
nf=0
n=0 # number of facilities, number of service areas
mthd=-1
if len(sys.argv) >= 2:
    fn=sys.argv[1]
if len(sys.argv) >= 3:
    maxr=float(sys.argv[2])
if len(sys.argv) >= 4:
    maxcp=int(sys.argv[3])
if len(sys.argv) >= 5:
    nf=int(sys.argv[4])
if len(sys.argv) >= 6:
    psize=int(sys.argv[5])
if len(sys.argv) >= 7:
    maxloops=int(sys.argv[6])
if len(sys.argv) >= 8:
    mthd=int(sys.argv[7])

def save_sol_location_problem():
    fn2=fn+"_sol"+str(n)+"_"+str(int(d.biobjective))+".txt"
    f = open(fn2,"w")
    idx=0
    f.write("Idx\tRID\tCenter\n")
    for i in range(d.num_units):
        k=d.node_groups[i]
        uk=d.facilityCandidate[k]
        di=d.nodes[i][4]
        dk=d.nodes[uk][4]
        s=str(di)+"\t"+str(k)+"\t"+str(dk)+"\n"
        f.write(s)
    f.close()

d.location_problem=1
d.all_units_as_candadate_locations=0
d.fixed_cost_obj=01
d.spatial_contiguity=0
d.adaptive_number_of_facilities=1
if nf>0: d.adaptive_number_of_facilities=0
d.pop_dis_coeff=100000#10000#1000000 #100000 for FLP 1-5 for PDP
d.pop_deviation=0.0 #for PDP
d.max_loops_solution_not_improved=maxloops #for search termination

d.initial_solution_method=0
#read instance file(s)
d.readfile(fn) #read self-defined instances
t0=time.time()
d.initial_solution_method=mthd  #0,1,2
# d.initial_solution_method=0 or 1 for a city with low-density poputation 
# d.initial_solution_method=2 for a city with high-density poputation
# d.initial_solution_method=-1 automatically select a method
d.solution_similarity_limit=5.0 #max(10.0,100.0/n)
d.solver_message=0
d.max_num_facility=nf
d.mip_solver="gurobi"#"cplex"
d.cflpr_matheuristic(maxr,maxcp,psize,maxloops)

d.search_stat()
print "=========================Final results========================="
print "objective:",d.biobjective
print "facility cost",d.objective_fcost
print "transportation cost:",d.objective
print "srrvice overload", d.objective_overload
print "pool size",len(d.region_pool)
print "total time",time.time()-t0
print "facilities selected",[x for x in d.centersID if x>=0]
print "demand assignment:", d.node_groups
print "service area stat:"
for i in range(d.num_districts):
    if d.district_info[i][0]==0: continue
    print d.facilityCandidate[i],d.district_info[i], d.district_info[i][2]/d.district_info[i][1],
    print "continuality?",d.check_continuality_feasibility(d.node_groups,i)

print "solutions: obj, f_cost,t_cost, overload"
for x in d.all_solutions: 
    print x[3]/d.total_pop,x[0],x[3],x[4],x[5]
    #print sum(d.facilityCapacity[y] for y in range(d.num_districts) if d.centersID[y]>=0),
    print [y for y in x[1] if y>=0]
    covered=[0 for i in range(20)]
    for i in range(d.num_units):
        k=x[2][i]#node_groups[i]
        dis=int(d.nodedij[i][k]*2)
        if dis<20: covered[dis]+=d.nodes[i][3]
    print "covered (%):",
    for i in range(20):
        if covered[i]==0: continue
        print [i,int(covered[i]*10000/d.total_pop)/100.0],  
    print
save_sol_location_problem()
print "-----------redius (km) coverage (%) statistics---------------"
covered=[0 for x in range(20)]
for i in range(d.num_units):
    k=d.all_solutions[0][2][i]#node_groups[i]
    dis=int(d.nodedij[i][k]*2)
    if dis<20: covered[dis]+=d.nodes[i][3]
for i in range(20):
    if covered[i]==0: continue
    print float(i+1)/2,"km", int(covered[i]*1000000/d.total_pop)/10000.0, "%"
print "---------------------------------------------------------------"
print

