# adjustable parameters
n <- 10 # total number of passengers
route_length <- n * 2 # origins and destinations
population_size <- 15 # population size for evolution
mutation_rate <- 1/n
termination_probe_length <- 200
debug_print_hillclimbing <- F

init <- function() {
    init_duration_matrix()
    init_population()
    evolution_history <<- c()

    cat("initial sequences:\n")
    print(population$seq)
    cat("initial fitness scores:\n")
    print(population$fitness)
}

evolve <- function(max_generations) {
    # approach: genetic algoritm, stop if not improving anymore
    gen <- 0
    while (gen < max_generations) {
        gen <- gen + 1
        next_gen()

        # progess check
        if (gen > termination_probe_length) {
            probe <- tail(evolution_history, termination_probe_length)
            older_half <- head(probe, termination_probe_length / 2)
            newer_half <- tail(probe, termination_probe_length / 2)
            if (sum(older_half) >= sum(newer_half)) {
                print("no progress, stopping")
                return()
            }
        }
    }
}

next_gen <- function() {
    # select parents
    parents <- select_parents()
    parent_sequences <- population$seq[parents,]

    # recombine and mutate
    sequence_is_duplicate <- T
    while (sequence_is_duplicate) {
        child_sequence <- recombine(parent_sequences)
        if (runif(1) < mutation_rate) {
            child_sequence <- mutate(child_sequence)
        }
        sequence_is_duplicate <- !is_new_sequence(child_sequence)
    }

    # form new population by replacing the least fit individual with the child
    child_route <- find_optimal_route(child_sequence)
    child_fitness <- target_function(child_route)
    swap_index <- which.min(population$fitness)
    update_individual(swap_index, child_sequence, child_route, child_fitness)

    # track progress
    mean_fitness <- mean(population$fitness)
    evolution_history <<- rbind(evolution_history, mean_fitness)
    generation_number <- length(evolution_history)
    cat("gen", generation_number, ": recombined", parents[1], parents[2], 
        "and put child at", swap_index, "-> mean fitness", mean_fitness, "\n")
}

best <- function() {
    fittest_index <- which.max(population$fitness)
    route <- population$route[fittest_index,]
    fitness <- population$fitness[fittest_index]
    cat("best route (", 1 / fitness, " minutes ): ", route)
}

# ----- INTERNAL FUNCTIONS ------

init_duration_matrix <- function() {
    # fill duration matrix with random duration values between 1 and 15
    duration_matrix <<- matrix(0, route_length, route_length)
    point_pairs <- combn(1:route_length, 2)
    for (i in 1:ncol(point_pairs)) {
        p1 <- point_pairs[1,i]
        p2 <- point_pairs[2,i]
        random_distance <- sample(15, 1)
        duration_matrix[p1,p2] <<- random_distance
        duration_matrix[p2,p1] <<- random_distance
    }
}

to_matrix_index <- function(p) if (p < 0) return(p + (2 * n) + 1) else return(p)
dur <- function(p1, p2) duration_matrix[to_matrix_index(p1), to_matrix_index(p2)]
orig <- function(f) f
dest <- function(f) -f

target_function <- function(route) {
    total_duration <- 0
    for (i in 2:length(route)) {
        start <- route[i-1]
        end <- route[i]
        total_duration <- total_duration + dur(start, end)
    }
    return(1 / total_duration) # higher values should indicate better solutions
}

# ----- outer optimization: best pickup sequence -----

sequence_slot_pairs <- combn(1:n, 2) # all permutations of length 2
sequence_slot_pairs_count <- ncol(sequence_slot_pairs)

init_population <- function() {
    population <<- list(
        seq=matrix(0,population_size,n), # pickup sequence
        route=matrix(,population_size,route_length), # calculated route for sequence
        fitness=vector(,population_size)) # calculated fitness for route
    
    # fill population with random sequences
    for (i in 1:population_size) {
        sequence_is_duplicate <- T
        while (sequence_is_duplicate) {
            sequence <- sample(n)
            sequence_is_duplicate <- !is_new_sequence(sequence)
        }
        route <- find_optimal_route(sequence)
        fitness <- target_function(route)
        update_individual(i, sequence, route, fitness)
    }
}

update_individual <- function(i, sequence, route, fitness) {
    population$seq[i,] <<- sequence
    population$route[i,] <<- route
    population$fitness[i] <<- fitness
}

is_new_sequence <- function(sequence) {
    matches <- apply(population$seq, 1, function(s) all(s == sequence))
    return(!any(matches))
}

select_parents <- function() {
    # approach: select random pair of parents, with probabilities proportional to fitness rank
    rand <- runif(1)
    rank <- 1
    while (rand < 0.66 && rank < (n-1)) {
        rank <- rank + 1
        rand <- rand * 2
    }
    fitness_ranks <- order(population$fitness, decreasing=T)
    selected_indices <- fitness_ranks[rank:(rank+1)]
    return(selected_indices)
}

recombine <- function(parents) {
    # approach: order crossover
    random_index <- sample(sequence_slot_pairs_count,1)
    random_pair <- sequence_slot_pairs[,random_index]
    slot_range <- random_pair[1]:random_pair[2]

    values_in_range <- parents[1,slot_range]
    fill_up_index <- 1
    recombined_sequence <- vector(,n)
    for (i in 1:n) {
        if (any(slot_range==i)) {
            recombined_sequence[i] <- parents[1,i]
        } else {
            value <- parents[2,fill_up_index]
            while (any(values_in_range==value)) {
                fill_up_index <- fill_up_index + 1
                value <- parents[2,fill_up_index]
            }
            recombined_sequence[i] <- value
            fill_up_index <- fill_up_index + 1
        }
    }
    return(recombined_sequence)
}

mutate <- function(sequence) {
    # approach: invert random partial sequence
    random_index <- sample(sequence_slot_pairs_count,1)
    start <- sequence_slot_pairs[1,random_index]
    end <- sequence_slot_pairs[2,random_index]
    mutated_sequence <- sequence
    mutated_sequence[start:end] <- rev(sequence[start:end])
    return(mutated_sequence)
}

# ----- inner optimization: best route for given pickup sequence -----

find_optimal_route <- function(sequence) {
    # approach: first calculate greedy solution, then optimize it locally
    greedy_solution <- greedy_route(sequence)
    optimized_solution <- optimize_route_locally(greedy_solution)
    return(optimized_solution)
}

greedy_route <- function(sequence) {
    # approach: insert dropoffs in the order of pickup
    route <- orig(sequence)

    for (i in 1:n) {
        passenger <- sequence[i]
        last_index <- length(route)
        pickup_index <- match(orig(passenger), route)
        if (pickup_index == last_index) {
            route <- append(route, dest(passenger), pickup_index)
            next
        }
        insert_positions <- pickup_index:last_index
        extended_routes <- t(sapply(insert_positions, function(pos) append(route, dest(passenger), pos)))
        target_values <- apply(extended_routes, 1, target_function)
        max_index <- which.max(target_values)
        route <- extended_routes[max_index,]
    }

    return(route)
}

optimize_route_locally <- function(route) {
    # approach: hillclimbing by moving one dropoff at a time within the route
    start_taget_value <- target_function(route)
    opt_route <- route
    opt_taget_value <- start_taget_value

    improved <- T
    rounds <- 0
    while (improved) {
        improved <- F
        rounds <- rounds + 1

        for (slot in 1:route_length) {
            stop <- route[slot]
            if (stop > 0) {
                next # stop is pickup
            }

            pickup_slot <- match(-stop, route)
            route_without_stop <- route[route != stop]

            for (insert_slot in pickup_slot:(route_length-1)) {
                variation <- append(route_without_stop, stop, insert_slot)

                target_value <- target_function(variation)
                if (target_value > opt_taget_value) {
                    opt_route <- variation
                    opt_taget_value <- target_value
                    improved <- T
                }
            }
        }
    }

    if (debug_print_hillclimbing && rounds > 1) {
        cat("optimized locally", rounds - 1, "times \n")
    }
    
    return(opt_route)
}