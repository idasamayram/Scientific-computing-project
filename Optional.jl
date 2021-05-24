### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 81440de5-1f02-4ae2-be61-9bb91c2c3c0e
begin
	import Pkg
	Pkg.add("VoronoiFVM")
	using VoronoiFVM
	Pkg.add("DifferentialEquations")
	using DifferentialEquations
	Pkg.add("LinearAlgebra")
	using LinearAlgebra
	Pkg.add("PyPlot")
	using Printf
	using PyPlot
end

# ╔═╡ 96eece23-a603-4c11-ba1f-65beb22a4fbb

import Base:push!


# ╔═╡ 882097cd-52e6-4733-b054-637ae2f09f21
md"""
# Optional Part
"""

# ╔═╡ 285f238a-9503-4f6c-9f82-1a85a72af321
md"""
## 1d Grid
"""

# ╔═╡ 9806327e-cd6d-4483-b285-c4b1767db35d
function barenblatt(x,t,m)
    tx=t^(-1.0/(m+1.0))
    xx=x*tx
    xx=xx*xx
    xx=1- xx*(m-1)/(2.0*m*(m+1));
    if xx<0.0
        xx=0.0
    end
    return tx*xx^(1.0/(m-1.0))
end

# ╔═╡ 53857250-a95f-458a-bf04-ed5b5fd32298
function create_porous_medium_problem(n,m,unknown_storage)
    h=1.0/convert(Float64,n/2)
    X=collect(-1:h:1)
    grid=VoronoiFVM.Grid(X)

    function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
    end

    storage!(f,u,node)= f[1]=u[1]

    physics=VoronoiFVM.Physics(flux=flux!,storage=storage!)

    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
    enable_species!(sys,1,[1])
    sys,X
end

# ╔═╡ 15702ae8-af6c-4d4c-981a-a21c11646843
function run_vfvm(;n=20,m=2,t0=0.001, tend=0.01,tstep=1.0e-6,unknown_storage=:dense)

    sys,X=create_porous_medium_problem(n,m,unknown_storage)

    inival=unknowns(sys)
    inival[1,:].=map(x->barenblatt(x,t0,m),X)

    solution=unknowns(sys)
    control=VoronoiFVM.NewtonControl()
    control.verbose=false
    control.Δt=tstep
    control.Δu_opt=0.05
    control.Δt_min=tstep

    times=collect(t0:t0:tend)
    times=[t0,tend]
    sol=VoronoiFVM.solve(inival,sys,times,control=control,store_all=true)
    err=norm(sol[1,:,end]-map(x->barenblatt(x,tend,m),X))
    sol,X,err
end

# ╔═╡ 0750c38b-4a6b-41ff-b926-ffd0972ac556
function run_diffeq(;n=20,m=2, t0=0.001,tend=0.01, unknown_storage=:dense,solver=nothing)
    sys,X=create_porous_medium_problem(n,m,unknown_storage)
    inival=unknowns(sys)
    inival[1,:].=map(x->barenblatt(x,t0,m),X)
    tspan = (t0,tend)
    sol=VoronoiFVM.solve(DifferentialEquations,inival,sys,tspan,solver=solver)
    err=norm(sol[1,:,end]-map(x->barenblatt(x,tend,m),X))
    sol, X,err
end

# ╔═╡ ae28420b-6917-40cb-9843-b419578ebf7f
function main(;m=2,n=20, solver=nothing, unknown_storage=:dense)
    function plotsol(sol,X)
        f=sol[1,:,:]'
        contourf(X,sol.t,f,0:0.1:10,cmap=:summer)
        contour(X,sol.t,f,0:1:10,colors=:black)
    end

    clf()
    subplot(121)

    t=@elapsed begin
        sol,X,err=run_vfvm(m=m,n=n, unknown_storage=unknown_storage)
    end
    title(@sprintf("VoronoiFVM: %.0f ms e=%.2e",t*1000,err))
    plotsol(sol,X)

    subplot(122)
    t=@elapsed begin
        sol,X,err=run_diffeq(m=m,n=n,solver=solver, unknown_storage=unknown_storage)
    end
    plotsol(sol,X)
    title(@sprintf("DifferentialEq: %.0f ms, e=%.2e",t*1000,err))

    gcf().set_size_inches(8,4)
    gcf()
end

# ╔═╡ aa870b11-1e32-44f9-be38-89d73a8503ae
main(;m=2,n=20, solver=nothing, unknown_storage=:dense)

# ╔═╡ 5d3b818f-6f98-4b2d-811a-362fba66bc5e
md"""
## 2d Grid
"""

# ╔═╡ d24bce0f-187d-4da8-b68c-f248f8157b7b
function create_porous_medium_problem2(n,m,unknown_storage)
    h=1.0/convert(Float64,n/2)
    X=collect(-1:h:1)
	Y=collect(-1:h:1)
    grid=VoronoiFVM.Grid(X,Y)

    function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
    end

    storage!(f,u,node)= f[1]=u[1]

    physics=VoronoiFVM.Physics(flux=flux!,storage=storage!)

    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
    enable_species!(sys,1,[1])
    sys,X
end

# ╔═╡ 50b8320c-a2f4-43f7-82b0-ccead70ef08c
function run_vfvm2(;n=20,m=2,t0=0.001, tend=0.01,tstep=1.0e-6,unknown_storage=:dense)

    sys,X,Y=create_porous_medium_problem(n,m,unknown_storage)

    inival=unknowns(sys)
    inival[1,:].=map((x,y)->barenblatt(x,t0,m),grid)

    solution=unknowns(sys)
    control=VoronoiFVM.NewtonControl()
    control.verbose=false
    control.Δt=tstep
    control.Δu_opt=0.05
    control.Δt_min=tstep

    times=collect(t0:t0:tend)
    times=[t0,tend]
    sol=VoronoiFVM.solve(inival,sys,times,control=control,store_all=true)
	err=norm(sol[1,:,end]-map((x,y)->barenblatt(x,tend,m),grid))
    sol,X,err
end

# ╔═╡ cc5c2c8d-24ae-445e-9fdf-090cd5b5f282
function run_diffeq2(;n=20,m=2, t0=0.001,tend=0.01, unknown_storage=:dense,solver=nothing)
    sys,X=create_porous_medium_problem(n,m,unknown_storage)
    inival=unknowns(sys)
    inival[1,:].=map((x,y)->barenblatt(x,t0,m),grid)
    tspan = (t0,tend)
    sol=VoronoiFVM.solve(DifferentialEquations,inival,sys,tspan,solver=solver)
	err=norm(sol[1,:,end]-map((x,y)->barenblatt(x,tend,m),grid))
    sol, X,err
end

# ╔═╡ 07a22413-86b7-4702-91c9-a82c4a167239
function main2(;m=2,n=20, solver=nothing, unknown_storage=:dense)
    function plotsol(sol,grid)
        f=sol[1,:,:]'
        contourf(grid,sol.t,f,0:0.1:10,cmap=:summer)
        contour(grid,sol.t,f,0:1:10,colors=:black)
    end

    clf()
    subplot(121)

    t=@elapsed begin
        sol,grid,err=run_vfvm(m=m,n=n, unknown_storage=unknown_storage)
    end
    title(@sprintf("VoronoiFVM: %.0f ms e=%.2e",t*1000,err))
    plotsol(sol,grid)

    subplot(122)
    t=@elapsed begin
        sol,grid,err=run_diffeq(m=m,n=n,solver=solver, unknown_storage=unknown_storage)
    end
    plotsol(sol,grid)
    title(@sprintf("DifferentialEq: %.0f ms, e=%.2e",t*1000,err))

    gcf().set_size_inches(8,4)
    gcf()
end

# ╔═╡ c28eb55f-a97d-42f7-a8b7-2bbcdbedbde2
main2(;m=2,n=20, solver=nothing, unknown_storage=:dense)

# ╔═╡ Cell order:
# ╠═81440de5-1f02-4ae2-be61-9bb91c2c3c0e
# ╠═96eece23-a603-4c11-ba1f-65beb22a4fbb
# ╟─882097cd-52e6-4733-b054-637ae2f09f21
# ╟─285f238a-9503-4f6c-9f82-1a85a72af321
# ╠═9806327e-cd6d-4483-b285-c4b1767db35d
# ╠═53857250-a95f-458a-bf04-ed5b5fd32298
# ╠═15702ae8-af6c-4d4c-981a-a21c11646843
# ╠═0750c38b-4a6b-41ff-b926-ffd0972ac556
# ╠═ae28420b-6917-40cb-9843-b419578ebf7f
# ╠═aa870b11-1e32-44f9-be38-89d73a8503ae
# ╟─5d3b818f-6f98-4b2d-811a-362fba66bc5e
# ╠═d24bce0f-187d-4da8-b68c-f248f8157b7b
# ╠═50b8320c-a2f4-43f7-82b0-ccead70ef08c
# ╠═cc5c2c8d-24ae-445e-9fdf-090cd5b5f282
# ╠═07a22413-86b7-4702-91c9-a82c4a167239
# ╠═c28eb55f-a97d-42f7-a8b7-2bbcdbedbde2
