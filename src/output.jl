function Plotres(N,pr,N_dt, output_folder :: String; save = false, name = "Nodal_pressure.png") #add more plots?
folder = joinpath(dirname(@__DIR__), "output", output_folder)
    if isdir(folder) == false
        mkdir(folder)
    end
fig_name =  joinpath(folder, string(output_folder,"_",name))
#clrs = [:red :blue :green  :cyan :magenta]
fig = plot(title="Nodal pressure", 
            xlabel="Time step", ylabel="Pressure [MPa]", 
            legend=:outerright, dpi = 300,
            size = (800,500))
    for n in N
        #plot!(1:N_dt, pr[:,n], lc=clrs[n], label="Node $n")
        plot!(1:N_dt, pr[:,n], label="Node $n")
    end
    display(fig)
    if save
        savefig(fig, fig_name) #how to set dpi higher
    end
end

function save_results(ES, output_folder :: String, name = "results.xlsx")
    folder = joinpath(dirname(@__DIR__), "output", output_folder)
    if isdir(folder) == false
        mkdir(folder)
    end
    xls_name = joinpath(folder, string(output_folder,"_",name))
    XLSX.openxlsx(xls_name, mode="w") do xf
    
        #write decision variables
        a = ES.out_dic
        # for k in keys(a)
        for k in eachindex(a)   
            sheet_name = string(a[k])            
            columns = value.(ES.model[a[k]])
            t = Tables.table(columns, header = collect(1:size(columns,2)) )
            if k == firstindex(a)
                sheet = xf[k]
                XLSX.rename!(sheet, string(a[k]))
                XLSX.writetable!(sheet, t)
                # XLSX.writetable(sheet, df)
            else
                sheet = XLSX.addsheet!(xf, string(a[k]))
                XLSX.writetable!(sheet, t)        
            end
        end
    
        #write feasibility_gap
        gap = calculate_feasibility_gap(ES)
        mae = mean(abs.(gap))
        mse = mean(gap.^2)
        rmse = sqrt.(mse)

        metrics = ["mae" mae; "mse" mse; "rmse" rmse]


        gap_t = Tables.table(gap, header = collect(1:size(gap,2)) )
        sheet = XLSX.addsheet!(xf, "feasibility_gap")
        XLSX.writetable!(sheet, gap_t)
    
        #anything else to write?
        # Parameters perhaps?
        # parsed_args
        sheet = XLSX.addsheet!(xf, "parameters")
        header = ["keys";"values"]
        #Convert symbols to strings
        args = Vector()
        for k in keys(parsed_args)
            push!(args, string(k))
        end
        
        par_data = [collect(args) collect(values(parsed_args))]
        par_data = [par_data; ["obj_val" objective_value(ES.model)]; ["solve_time" ES.solve_time] ; metrics]
        par_t = Tables.table(par_data, header = header)
        XLSX.writetable!(sheet, par_t)

    end
end