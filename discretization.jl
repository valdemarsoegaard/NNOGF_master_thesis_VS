"""
Created by Enrica Raheli and Yannick Werner
Licensed under xxxx
Language: Julia, version number
----------------------------------
Content of this file:
Function definitions for space discretization of a gas network.
"""

function space_discretization(nodes, pipes, segment)

    sub_nodes = nodes
    sub_pipelines = Dict{Int64, Any}()
    dx_vec=segment*ones(length(pipes))

    for (id, pipe) in pipes

        pipe_length = dx_vec[id] #segment length
        N_seg= div(pipe.LL, pipe_length) # number of segments of constant length
        if rem(pipe.LL, pipe_length) == 0 # remainder euqual to zero -> length is a multiple of discretization segment
            N_dx = N_seg # Number of segments in the pipe = number of segments of constant lenght
            last_length = pipe_length # length of the last subpipe
        else
            N_dx = N_seg+1 # Number of segments in the pipe = number of segments of constant lenght + one smaller segment
            last_length = rem(pipe.LL, pipe_length); # length of the last subpipe
        end

        # Calculate number of auxiliary subnodes between original nodes
        N_nodes_between = Int64(N_dx - 1)
    
        if N_nodes_between == 0
            #ToDo: Do we need that or can't we just use the original pipeline?
            pipeline_id = parse(Int64, string((pipe.start), (pipe.stop)))
            sub_pipelines[pipeline_id] = create_subpipe(
                pipe, "singular", pipeline_id, last_length, sub_nodes)
        else
            for n in 1:N_nodes_between
                ### Creation of (sub)nodes
                node_id = parse(Int64, string((pipe.start),(pipe.stop),n))
                sub_nodes[node_id] = create_subnode(nodes, pipe, node_id)
        
                ### Creation of (sub)pipelines
                if n == 1
                    pipeline_id = parse(Int64, string((pipe.start),(pipe.stop),n))
                    sub_pipelines[pipeline_id] = create_subpipe(
                        pipe, "first", pipeline_id, pipe_length, sub_nodes; segment_number=n, node_id=node_id)
                else
                    pipeline_id = parse(Int64, string((pipe.start),(pipe.stop),n))
                    sub_pipelines[pipeline_id] = create_subpipe(
                        pipe, "middle", pipeline_id, pipe_length, sub_nodes; segment_number=n, node_id=node_id)
                end
                
                if n == N_nodes_between
                    pipeline_id = parse(Int64, string((pipe.start),(pipe.stop),n+1))
                    sub_pipelines[pipeline_id] = create_subpipe(
                        pipe, "last", pipeline_id, last_length, sub_nodes; segment_number=n, node_id=node_id)
                end
                
            end
        end
    end
    return sub_nodes, sub_pipelines
end


function create_subnode(nodes, pipe, node_id)
    """
    Creates auxiliary subnodes between two original nodes 
    if the pipeline between them is discretized.
    """

    # The new subnodes don't have explicit pressure limits.
    # Hence, the minimum and maximum are used as dummy here.
    pmax_sub = max(nodes[pipe.start].pr_max,nodes[pipe.stop].pr_max)
    pmin_sub = min(nodes[pipe.start].pr_min,nodes[pipe.stop].pr_min)

    return Node(node_id, pmax_sub, pmin_sub, 0)
end


function create_subpipe(
    pipe,
    location,
    subpipeline_id,
    pipe_length,
    sub_nodes;
    segment_number=nothing,
    node_id=nothing)
    """
    Creates auxiliary subpipeline between two (sub-)nodes.
    """

    if location in ["first", "middle", "last"] == true
        if (isnothing(segment_number) || isnothing(node_id))
            throw(ArgumentError(
                "If location is 'first' or 'middle', segment_number must not be nothing"))
        end
    end

    # initialization of start and end node number
    start_node = nothing
    end_node = nothing

    if location == "singular"
        start_node = pipe.start
        end_node = pipe.stop
    elseif location in ["first", "middle"]
        start_node = ifelse(
            location == "first",
            pipe.start,
            parse(Int64, string(pipe.start, pipe.stop, segment_number-1)))
        end_node = node_id
    elseif location == "last"
        start_node = node_id
        end_node = pipe.stop
    else
        throw(ArgumentError(
            "Variable location must be one of 'first', 'middle', 'last', or 'singular."))
    end
    
    return Pipeline(subpipeline_id, start_node, end_node, pipe.fr, pipe.DD, pipe_length, pipe.cc, sub_nodes)

end