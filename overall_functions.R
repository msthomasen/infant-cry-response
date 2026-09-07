
plot_prior_posterior <- function(draws,
                                 posterior_name,
                                 prior_name,
                                 xlab,
                                 vline = NULL,
                                 title = NULL,
                                 central_prob = 0.99) {
  
  plot_df <- tibble(
    posterior = draws[, posterior_name],
    prior = draws[, prior_name]) %>%
    pivot_longer(
      everything(),
      names_to = "source",
      values_to = "value")
  
  limits <- quantile(
    plot_df$value,
    probs = c(
      (1 - central_prob) / 2,
      1 - (1 - central_prob) / 2),
    na.rm = TRUE)
  
  plot_df_trim <- plot_df %>%
    filter(
      value >= limits[1],
      value <= limits[2])
  
  ggplot(plot_df_trim, aes(value, fill = source, colour = source)) +
    geom_density(alpha = 0.25, linewidth = 0.5) +
    {
      if (!is.null(vline)) {
        geom_vline(xintercept = vline, linetype = 3)
      } else {
        NULL
      }} +
    theme_minimal() +
    labs(
      title = title,
      x = xlab,
      y = "Density",
      fill = NULL,
      colour = NULL)}


extract_summary <- function(samples, model_label, parameters) {
  samples$BUGSoutput$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    filter(parameter %in% parameters) %>%
    mutate(model = model_label) %>%
    dplyr::select(model, parameter, mean, sd, `2.5%`, `50%`, `97.5%`, Rhat)}

extract_theta_summary <- function(samples, model_label) {
  samples$BUGSoutput$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    filter(str_detect(parameter, "^theta\\[")) %>%
    mutate(
      from = str_match(parameter, "^theta\\[(\\d+),(\\d+)\\]")[, 2] %>% as.integer(),
      to = str_match(parameter, "^theta\\[(\\d+),(\\d+)\\]")[, 3] %>% as.integer()
    ) %>%
    transmute(
      model = model_label,
      from = factor(from, levels = 1:3, labels = c("none", "male", "female")),
      to = factor(to, levels = 1:3, labels = c("none", "male", "female")),
      mean = mean,
      lower = `2.5%`,
      upper = `97.5%`)}

extract_theta_draws <- function(samples, model_label) {
  as.matrix(as.mcmc(samples)) %>%
    as_tibble() %>%
    dplyr::select(matches("^theta\\[")) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      from = str_match(parameter, "^theta\\[(\\d+),(\\d+)\\]")[, 2] %>% as.integer(),
      to = str_match(parameter, "^theta\\[(\\d+),(\\d+)\\]")[, 3] %>% as.integer(),
      from = factor(from, levels = 1:3, labels = c("none", "male", "female")),
      to = factor(to, levels = 1:3, labels = c("none", "male", "female")),
      model = model_label)}



extract_theta_prior_draws <- function(samples, model_label) {
  as.matrix(as.mcmc(samples)) %>%
    as_tibble() %>%
    dplyr::select(matches("^theta_prior\\[")) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      from = str_match(parameter, "^theta_prior\\[(\\d+),(\\d+)\\]")[, 2] %>% as.integer(),
      to = str_match(parameter, "^theta_prior\\[(\\d+),(\\d+)\\]")[, 3] %>% as.integer(),
      from = factor(from, levels = 1:3, labels = c("none", "male", "female")),
      to = factor(to, levels = 1:3, labels = c("none", "male", "female")),
      model = model_label)}

extract_theta_p_prior_draws <- function(samples,
                                        model_label,
                                        dyads_keep = NULL,
                                        from_keep = NULL,
                                        to_keep = NULL) {
  draws_mat <- as.matrix(as.mcmc(samples))
  param_names <- colnames(draws_mat)
  
  theta_p_names <- param_names[str_detect(param_names, "^theta_p_prior\\[")]
  
  parsed <- tibble(parameter = theta_p_names) %>%
    mutate(
      dyad = str_match(parameter, "^theta_p_prior\\[(\\d+),(\\d+),(\\d+)\\]")[, 2] %>% as.integer(),
      from = str_match(parameter, "^theta_p_prior\\[(\\d+),(\\d+),(\\d+)\\]")[, 3] %>% as.integer(),
      to = str_match(parameter, "^theta_p_prior\\[(\\d+),(\\d+),(\\d+)\\]")[, 4] %>% as.integer())
  
  if (!is.null(dyads_keep)) {
    parsed <- parsed %>% filter(dyad %in% dyads_keep)}
  if (!is.null(from_keep)) {
    parsed <- parsed %>% filter(from %in% from_keep)}
  if (!is.null(to_keep)) {
    parsed <- parsed %>% filter(to %in% to_keep)}
  
  draws_mat %>%
    as_tibble() %>%
    dplyr::select(all_of(parsed$parameter)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    left_join(parsed, by = "parameter") %>%
    mutate(
      from = factor(from, levels = 1:3, labels = c("none", "male", "female")),
      to = factor(to, levels = 1:3, labels = c("none", "male", "female")),
      model = model_label)}


extract_theta_p_draws <- function(samples,
                                  model_label,
                                  dyads_keep = NULL,
                                  from_keep = NULL,
                                  to_keep = NULL) {
  draws_mat <- as.matrix(as.mcmc(samples))
  param_names <- colnames(draws_mat)
  
  theta_p_names <- param_names[str_detect(param_names, "^theta_p\\[")]
  
  parsed <- tibble(parameter = theta_p_names) %>%
    mutate(
      dyad = str_match(parameter, "^theta_p\\[(\\d+),(\\d+),(\\d+)\\]")[, 2] %>% as.integer(),
      from = str_match(parameter, "^theta_p\\[(\\d+),(\\d+),(\\d+)\\]")[, 3] %>% as.integer(),
      to = str_match(parameter, "^theta_p\\[(\\d+),(\\d+),(\\d+)\\]")[, 4] %>% as.integer())
  
  if (!is.null(dyads_keep)) {
    parsed <- parsed %>% filter(dyad %in% dyads_keep)}
  if (!is.null(from_keep)) {
    parsed <- parsed %>% filter(from %in% from_keep)}
  if (!is.null(to_keep)) {
    parsed <- parsed %>% filter(to %in% to_keep)}
  
  draws_mat %>%
    as_tibble() %>%
    dplyr::select(all_of(parsed$parameter)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    left_join(parsed, by = "parameter") %>%
    mutate(
      from = factor(from, levels = 1:3, labels = c("none", "male", "female")),
      to = factor(to, levels = 1:3, labels = c("none", "male", "female")),
      model = model_label)}

show_or_save_ggplot <- function(plot_png_path,
                                plot_object = NULL,
                                width = 8,
                                height = 5,
                                dpi = 300,
                                bg = "white") {
  if (file.exists(plot_png_path)) {
    img <- png::readPNG(plot_png_path)
    grid::grid.newpage()
    grid::grid.raster(img)
  } else {
    if (is.null(plot_object)) {
      stop("plot_object must be supplied when the PNG does not already exist.")}
    
    print(plot_object)
    
    ggsave(
      filename = plot_png_path,
      plot = plot_object,
      width = width,
      height = height,
      dpi = dpi,
      bg = bg)}}
