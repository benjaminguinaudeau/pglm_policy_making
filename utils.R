filter_only_if <- function(.tbl, criteria, ...){
  if(criteria){
    return(filter(.tbl, ...))
  } else {
    return(.tbl)
  }
}

get_electoral_cycle <- function(month, election_date, next_election, window){
  out <- case_when(
    as.numeric(month - lubridate::ymd(election_date)) <= window * 31 ~ "first_months",
    as.numeric(lubridate::ymd(next_election) - month) <= window * 31 ~ "last_months",
    T ~ "routine")
  return(out)
}  

sim_pglm <- function(fit, n_sim, new_data, form ){
  
  coef <- summary(fit)$estimate %>%
    as_tibble() %>%
    select_at(1:2) %>%
    set_names(c("coef", "sd"))
  
  k <- length(fit$estimate)
  
  V_beta <- vcov(fit) * array(coef$sd,c(k,k)) * t(array(coef$sd,c(k,k)))
  
  beta <- 1:n_sim %>%
    map_dfr(~{
      MASS::mvrnorm(1, coef$coef, V_beta) %>%
        imap_dfc(~tibble(.x) %>% set_names(.y))
    })
  
  req_cols <- str_split(as.character(form), "~|\\*|\\+") %>%
    reduce(c) %>%
    discard(~.x == "") %>%
    str_trim
  
  fit_data <- new_data %>%
    dplyr::select(req_cols) %>%
    mutate(.id = 1:nrow(.)) %>%
    drop_na
  
  sim_est <- new("sim",
                 coef = as.matrix(beta),
                 sigma = rep (sqrt(1), n_sim)) %>%
    coef %>%
    as_tibble %>%
    dplyr::select(colnames(modelr::model_matrix(fit_data, form))) %>%
    as.matrix()
  
  
  
  fit_data %>%
    modelr::model_matrix(form) %>%
    as.matrix %>%
    magrittr::multiply_by_matrix(t(sim_est)) %>%
    as_tibble %>%
    mutate_all(exp) %>%
    rename_all(~str_replace(.x, "V", "sim_")) %>%
    mutate(.id = fit_data$.id) %>%
    right_join(tibble(.id = 1:nrow(new_data))) %>%
    arrange(.id) %>%
    dplyr::select(-.id)
  
}

get_marginal_effect <- function(sims, term){
  
  terms <- str_split(term,  "\\*")[[1]] %>%
    str_trim
  
  tmp_dt <- sims %>%
    pivot_longer(cols = contains("sim_"), names_to = "sim_id", values_to = "expected_n_laws") %>%
    filter_at(vars(contains(terms[2])), ~!is.na(.x)) %>%
    get_quant_if_required(terms[2])
  
  tmp_dt %>%
    sample_n(10000) %>%
    ggplot(aes(x = !!dplyr::sym(terms[1]), y = expected_n_laws)) + 
    geom_jitter(alpha = .6) +
    geom_smooth(se = F, method = "lm") +
    theme_bw() +
    facet_wrap(as.formula(paste0("~",terms[2], "_quant")), ncol = 5) +
    ylim(c(0, 10))
}



get_quant_if_required <- function(tmp_dt, term){
  
  if(length(unique(tmp_dt[[term]])) > 5){
    tmp_dt <- tmp_dt %>%
      mutate(!!paste0(term, "_quant") := fct_reorder(paste0(term, "_", dvmisc::quant_groups(tmp_dt[[term]], 5)), tmp_dt[[term]]))
  } else {
    tmp_dt <- tmp_dt %>%
      mutate(!!paste0(term, "_quant") := paste0(term, "_", tmp_dt[[term]]))
    
  }
  return(tmp_dt)
}


get_marginal_effect_threeway <- function(sims, term){
  
  terms <- str_split(term,  "\\*")[[1]] %>% str_trim
  
  sal_term <- str_subset(terms, "salience")[1]
  other_term <- setdiff(terms, sal_term)
  
  sims %>%
    pivot_longer(cols = contains("sim_"), names_to = "sim_id", values_to = "expected_n_laws") %>%
    filter_at(vars(contains(other_term[1])), ~!is.na(.x)) %>%
    filter_at(vars(contains(other_term[2])), ~!is.na(.x)) %>%
    sample_n(10000) %>%
    get_quant_if_required(other_term[1]) %>%
    get_quant_if_required(other_term[2]) %>%
    dplyr::filter(expected_n_laws < 10) %>%
    ggplot(aes(x = !!dplyr::sym(sal_term), y = expected_n_laws)) + 
    geom_jitter(alpha = .6) +
    geom_smooth(se = F, method = "lm") +
    ylim(c(0, 10)) +
    theme_bw() +
    facet_grid(as.formula(paste0(paste0(other_term[1], "_quant"), "~",paste0(other_term[2], "_quant")))) 
  
}