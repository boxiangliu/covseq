library(data.table)
library(lubridate)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(countrycode)
library(scales)

metadata_fn = "../data/aggregated/metadata/merged.tsv"
out_dir = "../processed_data/metadata_statistics/"
dir.create(out_dir, recursive=TRUE)

read_metadata = function(metadata_fn){
    metadata = fread(metadata_fn)
    metadata = metadata[Collection_Date != "2022-04-05"]
    metadata = metadata[Collection_Date > "2019-12-01"]
    metadata = metadata[Country != ""]
    return(metadata)
}


process_metadata = function(metadata){
    parsed_date = parse_date_time(metadata[,Collection_Date], c("y", "ym", "ymd"))
    metadata$year_mon = sprintf("%s-%s",year(parsed_date), month(parsed_date))
    metadata$week = week(parsed_date)
    week_start = floor_date(parsed_date, unit="weeks")
    week_end = ceiling_date(parsed_date, unit="weeks")
    metadata$date_range = sprintf("%s to %s", week_start, week_end)
    metadata$week_start = week_start
    metadata$continent = countrycode(metadata$Country, origin="country.name", destination="continent")
    return(metadata)
}


get_sub_by_time = function(metadata){
    sub_by_time = metadata[, list(num=.N), by=c("week_start", "continent")]
    setorder(sub_by_time, 'week_start')
    sub_by_time[, cumsum := cumsum(num), by="continent"]
    sub_by_time_2 = sub_by_time[, list(continent="Total", num=sum(num)), by="week_start"]
    sub_by_time_2[, cumsum := cumsum(num)]
    sub_by_time = rbind(sub_by_time, sub_by_time_2)
    return(sub_by_time)
}


plot_sub_by_time = function(sub_by_time){
    p = ggplot(sub_by_time, aes(x=week_start,y=cumsum, color=continent)) + 
        geom_point() + 
        geom_line() + 
        scale_y_sqrt(labels=comma) + 
        scale_x_datetime(date_breaks="1 month", date_labels="%b") +
        scale_color_discrete(name = "Continent") + 
        ylab("Cumulative Number of Strains Collected") + 
        xlab("Time") + 
        theme(legend.position=c(0.05, 0.95), 
            legend.justification=c("left", "top"))
    return(p)
}


get_sub_by_country = function(metadata){
    sub_by_country = metadata[, list(submission=.N), by=c("Country", "continent")]
    setorder(sub_by_country, cols=-'submission')
    sub_by_country[,Country:=factor(Country, levels=Country)]
    return(sub_by_country)
}


plot_sub_by_country = function(sub_by_country){
    p = ggplot(sub_by_country, aes(x=Country, y=submission, fill=continent)) + 
        geom_bar(stat="identity") + 
        scale_y_log10(expand=c(0,0), name="Number of Submissions", labels=comma) + 
        scale_x_discrete(labels=NULL) + 
        scale_fill_discrete(name="Continent") +
        theme(axis.ticks.x = element_blank(), 
            legend.position = c(0.95, 0.95),
            legend.justification = c("right", "top"))
    return(p)
}


metadata = read_metadata(metadata_fn)
metadata = process_metadata(metadata)
sub_by_time = get_sub_by_time(metadata)
p1 = plot_sub_by_time(sub_by_time)
save_plot(sprintf("%s/sub_by_time.pdf", out_dir), p1, base_height=4, base_width=8)

sub_by_country = get_sub_by_country(metadata)
p2 = plot_sub_by_country(sub_by_country)
save_plot(sprintf("%s/sub_by_country.pdf", out_dir), p2, base_height=4, base_width=8)
