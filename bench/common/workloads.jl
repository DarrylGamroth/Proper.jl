module Workloads

export workload_matrix

function workload_matrix()
    return (
        grids=(256, 512, 1024),
        modes=(:python334, :corrected),
        backends=(:cpu, :cuda),
    )
end

end
