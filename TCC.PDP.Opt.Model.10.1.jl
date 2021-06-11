using JuMP, Gurobi, CSV, DataFrames

@time begin
    #-------------------------------------------------------tables
    tb_places = CSV.File("TimeDistance.csv") |> DataFrame
    tb_cars = CSV.File("Cars.csv") |> DataFrame
    tb_requests = CSV.File("Requests.csv") |> DataFrame
    
    #-------------------------------------------------------arrays
    requests = Array{Int64}(tb_requests[:,1])
    cars = Array{Int64}(tb_cars[:,1])
    nodes = Array{Int64}(1:max(maximum(tb_places[:,1]),maximum(tb_places[:,2])))
    
    #-------------------------------------------------------dics
    links = Tuple((tb_places[i,1],tb_places[i,2]) for i in 1:length(tb_places[:,1]))
    time_dict = Dict(links .=> tb_places[:,4] )
    distance_dict = Dict(links .=> tb_places[:,3])
    
    #-------------------------------------------------------model
    model = Model(with_optimizer(Gurobi.Optimizer, OutputFlag = 0))
    println("Graph.ok")
end

function PDP_BLD(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame, tb_cr::DataFrame, dist_dict::Dict)
        #-------------------------------------------------------variables
        @variables(m, begin
            w[rt], Bin
            x[nd,nd,cr], Bin
            y[nd,nd,cr,rt], Bin
            z[nd,nd,cr], Bin
            0.0 <= ta[nd,cr] <= 528.0
            0.0 <= td[nd,cr] <= 528.0
        end)
        
        #-------------------------------------------------------objective
        @objective(m, Max, sum(  tb_rt[r,11] * w[r] for r in rt) - 
                           sum(((dist_dict[i,j] * tb_cr[k,6]) * x[i,j,k] for i in nd for j in nd for k in cr)))
        println("Obj.ok")
end;

#-------------------------------------------------------constraint01
function PDP_C01(m::Model, nd::Array{Int64}, cr::Array{Int64}, tb_cr::DataFrame)
    for k in cr
        @constraint(m, ( sum( m[:x][tb_cr[k,2], j, k] for j in nd if j != tb_cr[k, 2] ) <= 1 ) ) 
    end
    println("c01.ok")
end;

#-------------------------------------------------------constraint02
function PDP_C02(m::Model, nd::Array{Int64}, cr::Array{Int64}, tb_cr::DataFrame)
    for k in cr
    @constraint(m, ( sum( m[:x][tb_cr[k,2], j, k] for j in nd if j != tb_cr[k,2] ) ==
                     sum( m[:x][j, tb_cr[k,3], k] for j in nd if j != tb_cr[k,3] ) 
                ) 
            ) 
    end                               
    println("c02.ok")
end;

#-------------------------------------------------------constraint03
function PDP_C03(m::Model, nd::Array{Int64}, cr::Array{Int64}, tb_cr::DataFrame)
    for i in nd, k in cr
        if i != tb_cr[k,2] && i != tb_cr[k,3] 
        @constraint(m, sum( m[:x][i,j,k] for j in nd if j != i ) - 
                       sum( m[:x][j,i,k] for j in nd if j != i ) == 0 )
        end   
    end
    println("c03.ok")
end;

#-------------------------------------------------------constraint04
function PDP_C04(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame)
    for r in rt 
        @constraint(m, sum(m[:y][tb_rt[r,2], j, k, r] for j in nd if j != tb_rt[r,2] for k in cr ) <= 1  )   
    end 
    println("c04.ok")
end;

#-------------------------------------------------------constraint05
function PDP_C05(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame)
    for r in rt
    @constraint(m, 
        (sum( m[:y][tb_rt[r,2],j, k, r] for j in nd if j != tb_rt[r,2] for k in cr) ==
         sum( m[:y][j,tb_rt[r,5], k, r] for j in nd if j != tb_rt[r,5] for k in cr))) 
    end                               
    println("c05.ok") 
end;

println("c06.out")

#-------------------------------------------------------constraint07
function PDP_C07(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame)
    for r in rt, i in nd
        if i != tb_rt[r,2] && i != tb_rt[r,5]
            for k in cr
                @constraint(m, sum( m[:y][i,j,k,r] for j in nd if j != i ) - 
                               sum( m[:y][j,i,k,r] for j in nd if j != i ) == 0 ) end
            end
    end
    println("c07.ok")
end;

#-------------------------------------------------------constraint08
function PDP_C08(m::Model, lks::Tuple, cr::Array{Int64}, rt::Array{Int64})
    for (i,j) in lks
        if j != i 
        for r in rt, k in cr
            @constraint(m, m[:y][i,j,k,r] <= m[:x][i,j,k] ) end
        end
    end
    println("c08.ok")
end;

#-------------------------------------------------------constraint09
function PDP_C09(m::Model, lks::Tuple, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame, tb_cr::DataFrame)
    for (i,j) in lks
        if j != i
            for k in cr
                @constraint(m,sum(tb_rt[r,8]*m[:y][ i, j, k, r] for r in rt) <= tb_cr[k,4] * m[:x][i,j,k]) end
            end
    end
        println("c09.ok")
end;

#-------------------------------------------------------constraint10
function PDP_C10(m::Model, lks::Tuple, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame, tb_cr::DataFrame)
    for (i,j) in lks
        if j != i
            for k in cr                                                     
                    @constraint(m,sum(tb_rt[r,9]*m[:y][ i, j, k, r] for r in rt) <= tb_cr[k,5] * m[:x][i,j,k]) end
            end
    end
    println("c10.ok")
end;

#-------------------------------------------------------constraint13
function PDP_C13(m::Model, lks::Tuple, cr::Array{Int64}, tb_cr::DataFrame)
    for (i,j) in lks
        if j != i
            for k in cr 
                if i != tb_cr[k,2] && j != tb_cr[k,3]
                    @constraint(m, m[:x][i,j,k] <= m[:z][i,j,k] ) end                                
            end
        end
    end
    println("c13.ok")
end;

#-------------------------------------------------------constraint14
function PDP_C14(m::Model, lks::Tuple, cr::Array{Int64}, tb_cr::DataFrame)
    for (i,j) in lks
        if j != i
            for k in cr 
                if i != tb_cr[k,2] && j != tb_cr[k,3]
                    @constraint(m, m[:z][i,j,k] + m[:z][j,i,k] == 1 ) end                                
            end
        end
    end
    println("c14.ok")
end;

#-------------------------------------------------------constraint15
function PDP_C15(m::Model, lks::Tuple, nd::Array{Int64}, cr::Array{Int64}, tb_cr::DataFrame)
    for (i,j) in lks
        if i != j
            for l in nd, k in cr 
                if j != l && l != i && i != tb_cr[k,2] && j != tb_cr[k,2] && l != tb_cr[k,3] 
                    @constraint(m, m[:z][i,j,k] + m[:z][j,l,k] + m[:z][l,i,k] <= 2 ) end    
            end
        end
    end
    println("c15.ok")
end;

#-------------------------------------------------------constraint16
function PDP_C16(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame)
    for r in rt 
        @constraint(m,sum(m[:y][i, tb_rt[r,5], k, r] for i in nd if i != tb_rt[r,5] for k in cr) == m[:w][r])   
    end 
    println("c16.ok")
end;

#-------------------------------------------------------constraint17
function PDP_C17(m::Model, lks::Tuple, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame)
    for (i,j) in lks
        if i != j
            for k in cr                                                      
                    @constraint(m,sum(tb_rt[r,10]*m[:y][ i, j, k, r] for r in rt) <= (1000000.0)) end
            end
        end
    println("c17.ok")
end;

#-------------------------------------------------------constraint18
function PDP_C18(m::Model, rt::Array{Int64})
    @constraint(m, sum( m[:w][r] for r in rt  ) >= (0.95) * size(rt,1) )   
    return println("c18.ok")
end;

#-------------------------------------------------------constraint19
function PDP_C19(m::Model, lks::Tuple, cr::Array{Int64}, t_dict::Dict )
    for (i,j) in lks
        if i != j
            for k in cr    
                    @constraint(m, ( ( m[:td][i,k] + t_dict[i,j] - m[:ta][j,k] ) <= ( (10000000000.0) * (1 - m[:x][i,j,k]) ) ) ) end
        end
    end
    println("c19.ok")
end; 

#-------------------------------------------------------constraint20
function PDP_C20(m::Model, nd::Array{Int64}, cr::Array{Int64})
    for k in cr 
        for i in nd
            @constraint(m, m[:ta][i,k] <= m[:td][i,k] )
        end
    end
    println("c20.ok")
end;

#-------------------------------------------------------constraint21
function PDP_C21(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame)
    for k in cars  
        for r in rt
            @constraints(m, begin
                            tb_rt[ r, 3 ] <= m[:ta][tb_rt[r,2], k ]
                            tb_rt[ r, 4 ] >= m[:td][tb_rt[r,2], k ]
            end)
        end
    end
    return println("c21.ok")
end;

#-------------------------------------------------------constraint22
function PDP_C22(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64}, tb_rt::DataFrame)
    for r in rt
        for k in cr
            @constraints(m, begin tb_rt[ r, 6 ] <= m[:ta][tb_rt[r,5], k ] 
                                  tb_rt[ r, 7 ] >= m[:td][tb_rt[r,5], k ] 
                end)
        end
    end
    return println("c22.ok")
end;

#-------------------------------------------------------constraint23
function PDP_C23(m::Model, cr::Array{Int64}, tb_cr::DataFrame)
    for k in cr
        @constraint(m, (m[:ta][tb_cr[k,3],k] - m[:td][tb_cr[k,2],k]) <= (528.0))
    end
    return println("c23.ok")
end;

#-------------------------------------------------------opt
function PDP_OPT(m::Model)
    JuMP.optimize!(m)
    println("opt.ok")
    return println("optimal solution = ", JuMP.objective_value(m))
end;

#-------------------------------------------------------obj
function PDP_SOL(m::Model)
    open("log.file.opt.txt", "w") do file write(file, "optimal solution = ", string(JuMP.objective_value(m))) end
    open("log.file.opt.txt", "a") do file write(file,"
", "------------------------------------------------","
")end
end; 

#-------------------------------------------------------w
function PDP_VRW(m::Model, rt::Array{Int64})
    for r in rt
        Rw = JuMP.value.(m[:w][r]) 
        if Rw == 1 (
        open("log.file.opt.txt", "a") do file write(file,"
", string((m[:w][r]))," = ",string(Rw)) end) end
    end
        open("log.file.opt.txt", "a") do file write(file,"
", "------------------------------------------------","
")end
end;

#-------------------------------------------------------x
function PDP_VRX(m::Model, nd::Array{Int64}, cr::Array{Int64})
    for k in cr
        for i in nd
            for j in nd
                Rx = JuMP.value.(m[:x][i,j,k])
                if Rx == 1 (
                open("log.file.opt.txt", "a") do file write(file,"
", string(m[:x][i,j,k])," = ",string(Rx)) end) end
            end
        end
    end          
open("log.file.opt.txt", "a") do file write(file,"
", "------------------------------------------------","
")end
end;

#-------------------------------------------------------y
function PDP_VRY(m::Model, nd::Array{Int64}, cr::Array{Int64}, rt::Array{Int64})
    for r in rt
        for k in cr
            for i in nd
                for j in nd 
                    Ry = JuMP.value.(m[:y][i,j,k,r])
                    if Ry == 1 (
                    open("log.file.opt.txt", "a") do file write(file,"
", string(m[:y][i,j,k,r])," = ",string(Ry)) end) end
                end
            end
        end
    end
open("log.file.opt.txt", "a") do file write(file,"
", "------------------------------------------------","
")end
end;

#-------------------------------------------------------z
function PDP_VRZ(m::Model, nd::Array{Int64}, cr::Array{Int64})
    for k in cr
        for i in nd
            for j in nd
                Rz = JuMP.value.(m[:z][i,j,k])
                if Rz == 1 (
                 open("log.file.opt.txt", "a") do file write(file,"
", string(m[:z][i,j,k])," = ",string(Rz))end) end
            end
        end
    end
open("log.file.opt.txt", "a") do file write(file,"
", "------------------------------------------------","
")end    
end;

#-------------------------------------------------------ta
function PDP_VTA(m::Model, nd::Array{Int64}, cr::Array{Int64})
        for k in cr
            for i in nd
                Rta = JuMP.value.(m[:ta][i,k])
                 open("log.file.opt.txt", "a") do file write(file,"
", string(m[:ta][i,k])," = ",string(Rta)) end
            end
        end
    open("log.file.opt.txt", "a") do file write(file,"
", "------------------------------------------------","
")   
    end
end;

#-------------------------------------------------------ta
function PDP_VTD(m::Model, nd::Array{Int64}, cr::Array{Int64})
    for k in cr
        for i in nd
            Rtd = JuMP.value.(m[:ta][i,k])
            open("log.file.opt.txt", "a") do file write(file,"
", string(m[:td][i,k])," = ",string(Rtd)) end
        end
    end
end;

BLD(model, nodes, cars, requests, tb_requests, tb_cars, distance_dict) =  
    @time PDP_BLD(model, nodes, cars, requests, tb_requests, tb_cars, distance_dict);
#---------------------------------------------------------
C01(model, nodes, cars, tb_cars)  =  @time PDP_C01(model, nodes, cars, tb_cars);
C02(model, nodes, cars, tb_cars)  =  @time PDP_C02(model, nodes, cars, tb_cars);
C03(model, nodes, cars, tb_cars)  =  @time PDP_C03(model, nodes, cars, tb_cars);
C04(model, nodes, cars, requests, tb_requests)  =  @time PDP_C04(model, nodes, cars, requests, tb_requests);
C05(model, nodes, cars, requests, tb_requests)  =  @time PDP_C05(model, nodes, cars, requests, tb_requests);
C07(model, nodes, cars, requests, tb_requests)  =  @time PDP_C07(model, nodes, cars, requests, tb_requests);
C08(model, links, cars, requests)  =  @time PDP_C08(model, links, cars, requests);
C09(model, links, cars, requests, tb_requests, tb_cars)  =  @time PDP_C09(model, links, cars, requests, tb_requests, tb_cars);
C10(model, links, cars, requests, tb_requests, tb_cars)  =  @time PDP_C10(model, links, cars, requests, tb_requests, tb_cars);
C13(model, links, cars, tb_cars)  =  @time PDP_C13(model, links, cars, tb_cars);
C14(model, links, cars, tb_cars)  =  @time PDP_C14(model, links, cars, tb_cars);
C15(model, links, nodes, cars, tb_cars)  =  @time PDP_C15(model, links, nodes, cars, tb_cars);
C16(model, nodes, cars, requests, tb_requests)  =  @time PDP_C16(model, nodes, cars, requests, tb_requests);
C17(model, links, cars, requests, tb_requests)  =  @time PDP_C17(model, links, cars, requests, tb_requests);
C18(model, requests)  =  @time PDP_C18(model, requests);
C19(model, links, cars, time_dict)  =  @time PDP_C19(model, links, cars, time_dict);
C20(model, nodes, cars)  =  @time PDP_C20(model, nodes, cars);
C21(model, nodes, cars, requests, tb_requests)  =  @time PDP_C21(model, nodes, cars, requests, tb_requests);
C22(model, nodes, cars, requests, tb_requests)  =  @time PDP_C22(model, nodes, cars, requests, tb_requests);
C23(model, cars, tb_cars)  =  @time PDP_C23(model, cars, tb_cars);
OPT(model)  =  @time PDP_OPT(model);
#---------------------------------------------------------
SOL(model)  =  @time PDP_SOL(model);
VRW(model, requests)  =  @time PDP_VRW(model, requests);
VRX(model, nodes, cars)  =  @time PDP_VRX(model, nodes, cars);
VRY(model, nodes, cars, requests)  =  @time PDP_VRY(model, nodes, cars, requests);
VRZ(model, nodes, cars) =  @time PDP_VRZ(model, nodes, cars);
VTA(model, nodes, cars)  =  @time PDP_VTA(model, nodes, cars);
VTD(model, nodes, cars)  =  @time PDP_VTD(model, nodes, cars);
print("fun.ok");

BLD(model, nodes, cars, requests, tb_requests, tb_cars, distance_dict);

C01(model, nodes, cars, tb_cars);

C02(model, nodes, cars, tb_cars);

C03(model, nodes, cars, tb_cars);

C04(model, nodes, cars, requests, tb_requests);

C05(model, nodes, cars, requests, tb_requests);

C07(model, nodes, cars, requests, tb_requests);

C08(model, links, cars, requests);

C09(model, links, cars, requests, tb_requests, tb_cars);

C10(model, links, cars, requests, tb_requests, tb_cars);

C13(model, links, cars, tb_cars);

C14(model, links, cars, tb_cars);

C15(model, links, nodes, cars, tb_cars);

C16(model, nodes, cars, requests, tb_requests);

C17(model, links, cars, requests, tb_requests);

C18(model, requests);

C19(model, links, cars, time_dict);

C20(model, nodes, cars);

C21(model, nodes, cars, requests, tb_requests);

C22(model, nodes, cars, requests, tb_requests);

C23(model, cars, tb_cars);

OPT(model);

SOL(model);

VRW(model, requests);

VRX(model, nodes, cars);

VRY(model, nodes, cars, requests);

VRZ(model, nodes, cars);

VTA(model, nodes, cars);

VTD(model, nodes, cars);
