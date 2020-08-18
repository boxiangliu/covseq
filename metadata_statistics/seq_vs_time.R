library(data.table)
library(lubridate)
library(cowplot)
library(countrycode)

metadata = fread("../data/aggregated/metadata/merged.tsv")
metadata = metadata[Collection_Date != "2022-04-05"]
metadata = metadata[Collection_Date > "2019-12-01"]
metadata = metadata[Country != ""]

parsed_date = parse_date_time(metadata[,Collection_Date], c("y", "ym", "ymd"))
metadata$year_mon = sprintf("%s-%s",year(parsed_date), month(parsed_date))
metadata$week = week(parsed_date)
week_start = floor_date(parsed_date, unit="weeks")
week_end = ceiling_date(parsed_date, unit="weeks")
metadata$date_range = sprintf("%s to %s", week_start, week_end)
metadata$week_start = week_start
metadata$continent = countrycode(metadata$Country, origin="country.name", destination="continent")

seq_vs_time = metadata[, list(num=.N), by=c("week_start", "continent")]
setorder(seq_vs_time, 'week_start')
seq_vs_time[, cumsum := cumsum(num), by="continent"]

ggplot(seq_vs_time, aes(x=week_start,y=cumsum, color=continent)) + 
  geom_point() + 
  geom_line() +
  ylab("Cumulative No. of Virus Collected") + 
  xlab("Time (Weekly Bars)") + 
  theme(axis.title=element_text(size=18,face="bold"))


sub_by_country = metadata[, list(submission=.N), by=c("Country", "continent")]
setorder(sub_by_country, cols=-'submission')
sub_by_country[,Country:=factor(Country, levels=Country)]

p = ggplot(sub_by_country, aes(x=Country, y=submission, fill=continent)) + 
    geom_bar(stat="identity") + 
    scale_y_log10() + 
    scale_x_discrete(labels=NULL) + 
    scale_fill_discrete(name="Continent") +
    theme(axis.ticks.x = element_blank(), 
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top")) + 
    annotate("text", x=2, y=55000, label="UK") + 
    annotate("text", x=4.5, y=35000, label="USA") + 
    annotate("text", x=6, y=9000, label="Australia") + 
    annotate("text", x=6, y=5000, label="India") + 
    annotate("text", x=7.5, y=3000, label="Spain")
print(p)

save_plot("test.pdf", p, base_height=4, base_width=8)
sub_by_country