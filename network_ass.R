library(igraph)
# --- Basic Network Generation Functions ---
generate.network.ER <- function(N, con) {
  A <- matrix(0, ncol = N, nrow = N)
  A[upper.tri(A)] <- rbinom(N * (N - 1) / 2, 1, con)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A
}

generate.network.B <- function(N, links.per.step) {
  L <- matrix(nrow = 0, ncol = 2)
  deg <- integer(N)
  for (i in 2:N) {
    n.new <- min(links.per.step, i - 1)
    linkto <- sample(i - 1, n.new, prob = deg[1:(i - 1)] + 1)
    newlinks <- cbind(rep(i, n.new), linkto)
    L <- rbind(L, newlinks)
    deg[i] = deg[i] + n.new
    deg[linkto] = deg[linkto] + 1
  }
  colnames(L) <- NULL
  L
}

# --- Core Epidemic Simulation Function ---
run_epidemic_and_monitor <- function(N, links, simlength, p.t, random_monitors, high_degree_monitors) {
  infected_status <- logical(N)
  patient_zero <- sample(N, 1)
  infected_status[patient_zero] <- TRUE
  
  infected_random_counts <- c(sum(infected_status[random_monitors]))
  infected_high_degree_counts <- c(sum(infected_status[high_degree_monitors]))
  
  for (t in 1:simlength) {
    discordant_links <- which(xor(infected_status[links[, 1]], infected_status[links[, 2]]))
    
    if (length(discordant_links) == 0) {
      infected_random_counts <- c(infected_random_counts, tail(infected_random_counts, 1))
      infected_high_degree_counts <- c(infected_high_degree_counts, tail(infected_high_degree_counts, 1))
      next
    }
    
    transmissions_occur <- rbinom(length(discordant_links), 1, p.t) == 1
    
    if (any(transmissions_occur)) {
      transmitting_edges <- links[discordant_links[transmissions_occur], ]
      nodes_involved <- unique(as.vector(transmitting_edges))
      newly_infected_nodes <- nodes_involved[!infected_status[nodes_involved]]
      infected_status[newly_infected_nodes] <- TRUE
    }
    
    infected_random_counts <- c(infected_random_counts, sum(infected_status[random_monitors]))
    infected_high_degree_counts <- c(infected_high_degree_counts, sum(infected_status[high_degree_monitors]))
  }
  return(list(random = infected_random_counts, high_degree = infected_high_degree_counts))
}

# --- Main Program: Setup and Run ---

# Parameters
N_nodes = 200
sim_duration = 30
trans_prob = 0.1
monitor_group_size = 20
num_epidemic_runs = 100

set.seed(Sys.time())

# 1. Generate Network
network_links <- generate.network.B(N_nodes, links.per.step = 2)
network_links
g_igraph <- graph_from_edgelist(network_links, directed = FALSE)
g_igraph

# 2. Identify Monitoring Groups
all_degrees <- degree(g_igraph)
sorted_by_degree <- order(all_degrees, decreasing = TRUE)
sorted_by_degree

random_monitors <- sample(1:N_nodes, monitor_group_size, replace = FALSE)
random_monitors
high_degree_monitors <- sorted_by_degree[1:monitor_group_size]
high_degree_monitors

# 3. Run Multiple Epidemic Simulations and Collect Data
all_runs_random_infected <- matrix(0, nrow = num_epidemic_runs, ncol = sim_duration + 1)
all_runs_high_degree_infected <- matrix(0, nrow = num_epidemic_runs, ncol = sim_duration + 1)

for (run in 1:num_epidemic_runs) {
  set.seed(run)
  run_results <- run_epidemic_and_monitor(N_nodes, network_links, sim_duration, trans_prob,
                                          random_monitors, high_degree_monitors)
  all_runs_random_infected[run, ] <- run_results$random
  all_runs_high_degree_infected[run, ] <- run_results$high_degree
}

#all_runs_random_infected
#all_runs_high_degree_infected

# 4. Calculate Averages
avg_random_infected <- colMeans(all_runs_random_infected)
avg_random_infected
avg_high_degree_infected <- colMeans(all_runs_high_degree_infected)
avg_high_degree_infected

# 5. Display Results
time_points <- 0:sim_duration
results_table <- data.frame(
  Time_Step = time_points,
  Avg_Infected_Random_Group = avg_random_infected,
  Avg_Infected_HighDegree_Group = avg_high_degree_infected
)

print(head(results_table, 5))
print(tail(results_table, 5))

# 6. Find "Early Warning" Time and Plot
threshold_infected = monitor_group_size * 0.2

time_random_reaches_threshold <- NA
for (i in 1:length(avg_random_infected)) {
  if (avg_random_infected[i] >= threshold_infected) {
    time_random_reaches_threshold <- time_points[i]
    break
  }
}

time_high_degree_reaches_threshold <- NA
for (i in 1:length(avg_high_degree_infected)) {
  if (avg_high_degree_infected[i] >= threshold_infected) {
    time_high_degree_reaches_threshold <- time_points[i]
    break
  }
}

plot(results_table$Time_Step, results_table$Avg_Infected_Random_Group, type = "l", col = "red", lwd = 2,
     xlab = "Time Step", ylab = "Avg Infected in Group",
     main = "Epidemic Spread: Random vs. High-Degree Monitors",
     ylim = c(0, monitor_group_size * 1.1))

lines(results_table$Time_Step, results_table$Avg_Infected_HighDegree_Group, col = "blue", lwd = 2)
legend("topleft", legend = c("Random Monitors", "High-Degree Monitors"),
       col = c("red", "blue"), lty = 1, lwd = 2)

abline(h = threshold_infected, lty = 2, col = "gray")

if (!is.na(time_random_reaches_threshold)) {
  points(time_random_reaches_threshold, threshold_infected, pch = 19, col = "red", cex = 1.5)
}
if (!is.na(time_high_degree_reaches_threshold)) {
  points(time_high_degree_reaches_threshold, threshold_infected, pch = 19, col = "blue", cex = 1.5)
}

if (!is.na(time_high_degree_reaches_threshold) && !is.na(time_random_reaches_threshold)) {
  early_warning_time <- time_random_reaches_threshold - time_high_degree_reaches_threshold
  text(mean(c(time_random_reaches_threshold, time_high_degree_reaches_threshold)), threshold_infected + 2,
       labels = paste("Early Warning:", early_warning_time, "steps"),
       cex = 0.8, col = "black")
}
time_random_reaches_threshold
time_high_degree_reaches_threshold
