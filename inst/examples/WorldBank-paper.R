library(animint2)
data(WorldBank)
WorldBank$region <- sub(" (all income levels)", "", WorldBank$region, fixed=TRUE)
wb.paper <- animint(
  title = "World Bank data viz for Animint paper",
  video = "https://vimeo.com/1087206179",
  source = "https://github.com/suhaani-agarwal/animint/blob/master/inst/examples/WorldBank-paper.R",
  ts = ggplot() +
    make_tallrect(WorldBank, "year") +
    guides(color = "none") +
    geom_line(
      aes(year, life.expectancy, group = country, colour = region),
      data = WorldBank,
      size = 4,
      alpha = 3/5,
      showSelected = "region",
      clickSelects = "country"
    ),
  scatter = ggplot() +
    geom_point(
      aes(fertility.rate, life.expectancy, colour = region, size = population, key = country),
      data = WorldBank,
      tooltip = paste(WorldBank$country, "population", WorldBank$population),
      key = "country",
      clickSelects = "country",
      showSelected = "year"
    ) +
    geom_text(
      aes(fertility.rate, life.expectancy, label = country),
      data = WorldBank,
      key = "country",
      clickSelects = "country",
      showSelected = c("country", "year", "region")
    ) +
    scale_size_animint(breaks = 10^(9:5)) +
    make_text(WorldBank, 5, 80, "year"),
  
  time = list(variable = "year", ms = 3000),
  duration = list(year = 1000),
  selector.types = list(
    country = "multiple",
    region = "multiple"
  ),
  first = list(
    year = 1979,
    country = c("United States", "Vietnam"),
    region = c("East Asia & Pacific", "North America")
  )
)

animint2pages(wb.paper, "WorldBank-paper")