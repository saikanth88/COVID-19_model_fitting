library(tidyverse)
library(lubridate)
library(dplyr)
library(viridis)
library(gridExtra)
library(grid)
library(ggsn)

PATH_TO_FILES = "C:/users/saika/desktop/Grid_Search_COVID/outputs/"

county = c("high", "mid", "low")

plot_pops = vector("list", length = length(county))
plot_beta = vector("list", length = length(county))
plot_hosp = vector("list", length = length(county))

#----------------------------------------------------------------------------------------------
plot_theme = list(
  labs(x = "", y = ""),
  theme_classic(),
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    title = element_text(size = 8, color = "black"),
    legend.position = "none",
    legend.direction = "horizontal",
    plot.margin = margin(t=5.5, r=5.5, b=12, l=5.5)
  )
)

region_col = c("#F05555", "#3CB371", "#6495ED")

#----------------------------------------------------------------------------------------------

for(i in 1:length(county)){
  
  this_county = county[i]
  
  pop_local_df = read_csv(paste0(PATH_TO_FILES, this_county, "_countyData.csv"))
  
  plot_pops[[i]] = ggplot() +
    scale_size(name = "Pop.Size", limits = c(200,3300000),  range = c(1, 15), 
               breaks = c(100, 1000, 10000, 100000, 1000000)) +
    geom_point(data = pop_local_df, aes(x = long, y = lat, size = pop_size), 
               alpha = 1.0, shape = 21, fill = region_col[i]) +
    coord_quickmap() +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          title = element_text(size = 7, color = "black"),
          legend.position = "none",
          plot.margin = margin(t=5.5, r=5.5, b=30, l=5.5)) +
    ggsn::scalebar(data = pop_local_df,
                   x.min = range(long)[1],
                   x.max = range(long)[2],
                   y.min = range(lat)[1],
                   y.max = range(lat)[2],
                   location = "bottomright",
                   anchor = c(x = range(pop_local_df$long)[2] - 0.02, y = range(pop_local_df$lat)[1] - 0.25),
                   dist = 40, dist_unit = "km",
                   transform = TRUE,
                   st.size = 2.5,
                   st.dist = 0.05,
                   height = 0.05,
                   border.size = 0.4)
  
  if(i == 2){ 
    plot_pops[[i]] = plot_pops[[i]] + 
      theme(legend.position = c(0.5,-0.6), legend.direction = "horizontal")}
  
  if(i == 3){ 
    plot_pops[[i]] = plot_pops[[i]] + 
      scale_x_continuous(limits = c(NA, range(pop_local_df$long)[2] + 0.2))}
  
  model_data_df = read_csv(paste0(PATH_TO_FILES, "beta_group/", this_county, "_params_all_counts_5000_all_pops.csv"))
  n_days = max(model_data_df$day) 
  
  start_date = mdy("03-01-2020")
  end_date = start_date + n_days
  
  date_seq = seq(start_date, end_date, by = "day")
  date_breaks = seq(start_date, end_date, by = "1 month")
  
  model_data_df = model_data_df %>%
    group_by(iter, realz) %>%
    mutate(date = date_seq) %>%
    ungroup() 
  
  raw_data_df =  read_csv(paste0(PATH_TO_FILES, this_county, "_R0_hosp_count.csv"))
  
  raw_data_df = raw_data_df %>%
    group_by(day) %>%
    summarize(n_hosp = sum(n_hosp))
  
  raw_data_df = raw_data_df %>%
    filter(day <= (n_days+1)) %>%
    mutate(date = date_seq)
  
  grid_param_df =  read_csv(paste0(PATH_TO_FILES, "beta_group/", this_county, "_full_covid_list_5000_all_pops.csv"), 
                            col_names = FALSE)
  
  n_param = 1
  n_col   = ncol(grid_param_df)
  n_weeks = (n_col - 1)/n_param
  
  colnames(grid_param_df) = c(sprintf("beta_%02d",seq(01:n_weeks)),
                              "Lhood")
  
  n_top = 200
  
  top_param_df = top_n(grid_param_df, n_top, grid_param_df[n_col])
  
  top_beta_df = top_param_df[, c(1:n_weeks)] %>%
    gather(key = "beta",value = "value", beta_01:beta_40) %>%
    group_by(beta) %>%
    summarize(low_beta  = quantile(value, 0.025),
              med_beta  = quantile(value, 0.5),
              high_beta = quantile(value, 0.975))

  qtl_beta_df = top_beta_df[rep(seq_len(nrow(top_beta_df)), each = 7),]%>%
    mutate(date = date_seq)

  raw_beta_df =  read_csv(paste0(PATH_TO_FILES, this_county, "_beta_vals.csv"))
  
  raw_beta_df = raw_beta_df %>%
    filter(day <= 7*n_weeks) %>%
    group_by(realz, day) %>%
    summarize(beta = mean(beta))
  
  raw_beta_df = raw_beta_df %>%
    mutate(date = date_seq) %>%
    mutate(week = cut(date, breaks = "7 days", labels = FALSE)) %>%
    group_by(week) %>%
    mutate(avg_beta = mean(beta))
  
  qtl_n_hosp = model_data_df %>%
    group_by(date) %>%
    summarize(low_n_hosp  = quantile(n_hosp, 0.025),
              med_n_hosp  = quantile(n_hosp, 0.5),
              high_n_hosp = quantile(n_hosp, 0.975))
  
  plot_beta[[i]] =
    ggplot() + plot_theme +
    labs(title = "Weekly transmission rate") +
    scale_x_date(breaks = date_breaks, date_labels = "%b-%Y") +
    geom_ribbon(data = qtl_beta_df, aes(x = date, ymin = low_beta, ymax = high_beta),
                fill = "#D16103", alpha = 0.5) +
    geom_path(data = qtl_beta_df, aes(x = date, y = med_beta, color = "model"), size = 0.75) +
    geom_path(data = raw_beta_df, aes(x = date, y = avg_beta, color = "data"), size = 0.5) +
    scale_color_manual(name = " ", values = c("model" = "#D16103", "data" = "#778899"))

  plot_hosp[[i]] = 
    ggplot() + plot_theme + 
    labs(title = "New hospitalizations per day") +
    scale_x_date(breaks = date_breaks, date_labels = "%b-%Y") +
    geom_ribbon(data = qtl_n_hosp, aes(x = date, ymin = low_n_hosp, ymax = high_n_hosp),
                fill = "#00AFBB", alpha = 0.5) +
    geom_path(data = qtl_n_hosp, aes(x = date, y = med_n_hosp, color = "model"), size = 1) +
    geom_path(data = raw_data_df, aes(x = date, y = n_hosp, color = "data"), size = 0.5) +
    scale_color_manual(name = " ", values = c("model" = "#00AFBB", "data" = "#778899"))
  
}

plot_beta[[2]] = plot_beta[[2]]+ theme(legend.position = c(0.5,-0.4))
plot_hosp[[2]] = plot_hosp[[2]]+ theme(legend.position = c(0.5,-0.4))

plot_combo_pop =
  arrangeGrob(arrangeGrob(plot_pops[[1]], nrow = 1,
                          top = textGrob("Urban",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              arrangeGrob(plot_pops[[2]], nrow = 1,
                          top = textGrob("Suburban",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              arrangeGrob(plot_pops[[3]], nrow = 1,
                          top = textGrob("Rural",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              ncol=3)

plot_combo_beta =
  arrangeGrob(arrangeGrob(plot_beta[[1]], nrow = 1,
                          top = textGrob("Urban",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              arrangeGrob(plot_beta[[2]], nrow = 1,
                          top = textGrob("Suburban",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              arrangeGrob(plot_beta[[3]], nrow = 1,
                          top = textGrob("Rural",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              ncol=3)

plot_combo_hosp =
  arrangeGrob(arrangeGrob(plot_hosp[[1]],nrow = 1,
                          top = textGrob("Urban",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              arrangeGrob(plot_hosp[[2]], nrow = 1,
                          top = textGrob("Suburban",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              arrangeGrob(plot_hosp[[3]], nrow = 1,
                          top = textGrob("Rural",
                                         gp = gpar(fontsize = 10, fontface = "bold"))),
              ncol=3)

ggsave(filename = paste0("regions", ".png"),
       plot = plot_combo_pop,
       path = paste0(PATH_TO_FILES, "beta_group/"),
       device = "png",
       height = 3,
       width = 8,
       units = "in",
       dpi = 300)

ggsave(filename = paste0("beta_estimation", ".png"),
       plot = plot_combo_beta,
       path = paste0(PATH_TO_FILES, "beta_group/"),
       device = "png",
       height = 3,
       width = 8,
       units = "in",
       dpi = 300)

ggsave(filename = paste0("hosp_fitting", ".png"),
       plot = plot_combo_hosp,
       path = paste0(PATH_TO_FILES, "beta_group/"),
       device = "png",
       height = 3,
       width = 8,
       units = "in",
       dpi = 300)