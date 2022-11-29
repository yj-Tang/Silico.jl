using Plots
using Statistics
using Random
using JLD2
using CUDA
using Flux
using BSON
using BenchmarkTools
CUDA.functional()

include("methods.jl")



function extract_feature_label(mechanism, storage::TraceStorage{T,H}, i) where {T,H}
    @assert i > 1
    solver = mechanism.solver
    data = solver.data
    problem = solver.problem
    indices = solver.indices
    methods = solver.methods
    cone_methods = solver.cone_methods
    solution = solver.solution

    idx_duals = indices.duals
    idx_slacks = indices.slacks
    idx_equality = indices.equality
    idx_slackness = indices.slackness

    previous_variables = storage.variables[i-1]
    variables = storage.variables[i]
    previous_parameters = storage.parameters[i-1]
    parameters = storage.parameters[i]
    solution.all .= previous_variables

    # compute the residual of the previous solution under the current parameters
    Mehrotra.evaluate!(problem, methods, cone_methods, solution, parameters,
        equality_constraint=true,
        equality_jacobian_variables=true,
        equality_jacobian_parameters=true,
        cone_constraint=true,
        cone_jacobian=false,
        compressed=false,
        sparse_solver=false)
    Mehrotra.residual!(data, problem, indices,
        residual=true,
        jacobian_variables=true,
        jacobian_parameters=true,
        compressed=false,
        sparse_solver=false)

    unoptimized_residual = deepcopy(data.residual.all)
    # previous_jacobian = vec(deepcopy(data.jacobian_variables_compressed_dense))
    # previous_log_variables = log.(10, previous_variables[[idx_duals; idx_slacks]])
    previous_active_set = previous_variables[idx_duals] .>= previous_variables[idx_slacks]
    active_set = variables[idx_duals] .>= variables[idx_slacks]
    # δ_parameters = parameters - previous_parameters
    # δ_residual = data.jacobian_parameters * δ_parameters

    # xi = [previous_active_set;
    #     # previous_log_variables;
    #     # previous_jacobian;
    #     previous_variables;
    #     parameters;
    #     unoptimized_residual[idx_equality]; # the idx_complementarity always = 0
    #     δ_residual[idx_equality]; # the idx_complementarity always = 0
    #     δ_parameters]


    solution.all .= 0.0
    Mehrotra.evaluate!(problem, methods, cone_methods, solution, parameters,
        equality_constraint=true,
        equality_jacobian_variables=true,
        equality_jacobian_parameters=true,
        cone_constraint=true,
        cone_jacobian=false,
        compressed=false,
        sparse_solver=false)
    Mehrotra.residual!(data, problem, indices,
        residual=true,
        jacobian_variables=true,
        jacobian_parameters=true,
        compressed=false,
        sparse_solver=false)
    zero_residual = deepcopy(data.residual.all)
    sdf = - zero_residual[idx_slackness][[1,5]]

    anticipated_sdf = - unoptimized_residual[idx_slackness][[1,5]] + previous_variables[idx_slacks][[1,5]]
    previous_pose = parameters[1:3]
    previous_velocity = parameters[4:6]
    input = parameters[7:9]
    xi = [
        # active_set[[1,5]]; # impact only
        previous_active_set[[1,5]]; # impact only
        previous_pose;
        previous_velocity;
        input;
        # active_set[[1,5]]; # impact only
        sdf; # sdf(current)
        anticipated_sdf; # sdf(current + Δt vel)
    ]
    yi = active_set[[1,5]]
    return xi, yi
end


mech.solver.solution.primals
mech.solver.solution.duals
mech.solver.parameters
mech.bodies[1]

################################################################################
# load data
################################################################################
batch_size = 500
x_train, y_train, x_val, y_val, x_test, y_test, μ, σ = load_dataset(; name="composed_capsule_dataset_0")
n_input = size(x_train, 1)
n_output = size(y_train, 1)

train_loader = Flux.DataLoader((x_train, y_train), batchsize=batch_size, shuffle=true)
x_train = gpu(x_train)
y_train = gpu(y_train)
x_val = gpu(x_val)
y_val = gpu(y_val)
x_test = gpu(x_test)
y_test = gpu(y_test)
μ = gpu(μ)
σ = gpu(σ)

################################################################################
# define models
################################################################################
baseline_model(x) = ((x .* (1e-5 .+ σ)) .+ μ)[1:n_output,:]
baseline_model(x_train)

cpu_model = Chain(
    Dense(n_input => 80, tanh),
    Dense(80 => 80, tanh),
    Dense(80 => 50, tanh),
    Dense(50 => 50, tanh),
    Dense(50 => 50, tanh),
    Dense(50 => 30, tanh),
    Dense(30 => n_output, sigmoid))
model = fmap(cu, cpu_model)
parameters = Flux.params(model)


################################################################################
# define losses
################################################################################
loss(x, y) = Flux.Losses.mse(model(x), y)
baseline_loss(x, y) = Flux.Losses.mse(baseline_model(x), y)

loss(x_train, y_train)
baseline_loss(x_train, y_train)


################################################################################
# training
################################################################################
n_epoch = 101
optimizer = Adam(0.001, (0.9, 0.999), 1.0e-8)
validation_loss() = round(loss(x_val, y_val), digits=4)

train_model!(train_loader, loss, parameters, optimizer, n_epoch;
    validation_loss=validation_loss,
    print_epoch=5)


save_model(model, name="composed_capsule_model_0")
loaded_cpu_model = load_model(name="composed_capsule_model_0")


################################################################################
# Analysis
################################################################################
loss(x_train, y_train)
loss(x_val, y_val)
loss(x_test, y_test)
baseline_loss(x_train, y_train)
baseline_loss(x_val, y_val)
baseline_loss(x_test, y_test)

sum(error_distribution(x_train, y_train, m=model)[1:1]) / size(x_train,2)
sum(error_distribution(x_val,   y_val,   m=model)[1:1]) / size(x_val,2)
sum(error_distribution(x_test,  y_test,  m=model)[1:1]) / size(x_test,2)

sum(error_distribution(x_train, y_train, m=model)[1:2]) / size(x_train,2)
sum(error_distribution(x_val,   y_val,   m=model)[1:2]) / size(x_val,2)
sum(error_distribution(x_test,  y_test,  m=model)[1:2]) / size(x_test,2)

sum(error_distribution(x_train, y_train, m=model)[1:3]) / size(x_train,2)
sum(error_distribution(x_val,   y_val,   m=model)[1:3]) / size(x_val,2)
sum(error_distribution(x_test,  y_test,  m=model)[1:3]) / size(x_test,2)

error_distribution(x_train, y_train, m=baseline_model)[1] / size(x_train,2)
error_distribution(x_val, y_val, m=baseline_model)[1] / size(x_val,2)
error_distribution(x_test, y_test, m=baseline_model)[1] / size(x_test,2)

error_distribution(x_test,  y_test, m=baseline_model)
error_distribution(x_test,  y_test, m=model)
binary_projection(x_train, y_train, 0.05, m=model)