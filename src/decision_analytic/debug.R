df <- data.frame(matrix(nrow=5, ncol = 2))

names(df) <- c("variable", "percentage")
df$variable <- c("Carbohydrates", "Warming", "NGTnotPresent", "DrainNotPresent", "DrEaMing")
df$percentage <- c(0.67,0.33,0.86,0.78,0.58)

df <- df %>% mutate(group=ifelse(percentage <0.6, "red",
                                 ifelse(percentage>=0.6 & percentage<0.8, "orange","green")),
                    label=paste0(percentage*100, "%"),
                    title=dplyr::recode(variable, `Carbohydrates`="Preoperative\ncarbohydrate loading",
                                        `Warming`="Intraoperative\nwarming",
                                        `NGTnotPresent`="Patients without a\nnasogastric tube\non arrival in recovery",
                                        `DrainNotPresent`="Patients without an\nabdominal drain\non arrival in recovery",
                                        `DrEaMing`="Patients DrEaMing on\npostoperative day 1"))


ggplot(df, aes(fill = group, ymax = percentage, ymin = 0, xmax = 2, xmin = 1)) +
  geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="#ece8bd") +
  geom_rect() + 
  coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2)) + ylim(c(0,2)) +
  geom_text(aes(x = 0, y = 0, label = label, colour=group), size=6.5, family="Poppins SemiBold") +
  geom_text(aes(x=1.5, y=1.5, label=title), family="Poppins Light", size=4.2) + 
  facet_wrap(~title, ncol = 5) +
  theme_void() +
  scale_fill_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188")) +
  scale_colour_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill=FALSE) +
  guides(colour=FALSE)
