install.packages("treemap")
library(treemap)

data(GNI2014)
head(GNI2014)

treemap(GNI2014,
        index="continent", #指定分组的列
        vSize="population", #指定面积大小的列
        vColor="GNI", #指定颜色深浅的列
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = FALSE, big.mark = " "))


# colors indicate density (like a population density map)
treemap(GNI2014,
        index=c("continent","country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="population", #指定面积大小的列
        vColor="GNI", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
install.packages("readxl")
library(readxl)


# Officially my code. Can run from here

# Officially my code. Can run from here
library(treemap)
library(readxl)
data_pop_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_pop_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_pop_neg_2100.pdf", width = 4, height = 5)
treemap(data_pop_neg,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Population", #指定面积大小的列
        vColor="Negatively-Impacted Population", #指定颜色深浅的列
        type="dens")
        #title.legend = NA,  # Remove legend title
        #position.legend = "none")  # Turn off legend)
dev.off()

data_pop_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_pop_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_pop_neg_2100_Hunger.pdf", width = 4, height = 5)
treemap(data_pop_neg,
        index=c("Hunger","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Population", #指定面积大小的列
        vColor="Negatively-Impacted Population", #指定颜色深浅的列
        type="dens")
        #title.legend = NA,  # Remove legend title
        #position.legend = "none")  # Turn off legend)
dev.off()


data_pop_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_pop_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_pop_neg_2100_FSI.pdf", width = 4, height = 5)
treemap(data_pop_neg,
        index=c("FSI","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Population", #指定面积大小的列
        vColor="Negatively-Impacted Population", #指定颜色深浅的列
        type="dens")
        #title.legend = NA,  # Remove legend title
        #position.legend = "none")  # Turn off legend)
dev.off()


data_pop_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_pop_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_pop_neg_2100_Gender.pdf", width = 4, height = 5)
treemap(data_pop_neg,
        index=c("Gender","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Population", #指定面积大小的列
        vColor="Negatively-Impacted Population", #指定颜色深浅的列
        type="dens")
        #title.legend = NA,  # Remove legend title
        #position.legend = "none")  # Turn off legend)
dev.off()


data_pop_posi <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_pop_posi.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_pop_posi_2100.pdf", width = 4, height = 4)
treemap(data_pop_posi,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Positively-Impacted Population", #指定面积大小的列
        vColor="Positively-Impacted Population", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_cattle_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_cattle_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_cattle_neg_2100.pdf", width = 4, height = 4)
treemap(data_cattle_neg,
        index = c("Development", "Country"),  # Specify grouping
        vSize = "Negatively-Impacted Cattle", # Specify size
        vColor = "Negatively-Impacted Cattle", # Specify color
        type = "dens") 
        #title.legend = NA,  # Remove legend title
        #position.legend = "none")  # Turn off legend
dev.off()

data_cattle_posi <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_cattle_posi.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_cattle_posi_2100.pdf", width = 4, height = 4)
treemap(data_cattle_posi,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Positively-Impacted Cattle", #指定面积大小的列
        vColor="Positively-Impacted Cattle", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_sheep_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_sheep_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_sheep_neg_2100.pdf", width = 4, height = 4)
treemap(data_sheep_neg,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Sheep", #指定面积大小的列
        vColor="Negatively-Impacted Sheep", #指定颜色深浅的列
        type="dens")
        #title.legend = NA,  # Remove legend title
        #position.legend = "none")  # Turn off legend)
dev.off()

data_sheep_posi <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_sheep_posi.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_sheep_posi_2100.pdf", width = 4, height = 4)
treemap(data_sheep_posi,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Positively-Impacted Sheep", #指定面积大小的列
        vColor="Positively-Impacted Sheep", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_goats_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_goats_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_goats_neg_2100.pdf", width = 4, height = 4)
treemap(data_goats_neg,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Goats", #指定面积大小的列
        vColor="Negatively-Impacted Goats", #指定颜色深浅的列
        type="dens")
       # title.legend = NA,  # Remove legend title
       # position.legend = "none")  # Turn off legend)
dev.off()

data_goats_posi <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_goats_posi.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_goats_posi_2100.pdf", width = 4, height = 4)
treemap(data_goats_posi,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Positively-Impacted Goats", #指定面积大小的列
        vColor="Positively-Impacted Goats", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_livestock_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_livestock_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_livestock_neg_2100.pdf", width = 4, height = 4)
treemap(data_livestock_neg,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Livestock", #指定面积大小的列
        vColor="Negatively-Impacted Livestock", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_livestock_posi <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_livestock_posi.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_livestock_posi_2100.pdf", width = 4, height = 4)
treemap(data_livestock_posi,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Positively-Impacted Livestock", #指定面积大小的列
        vColor="Positively-Impacted Livestock", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()


library(treemap)
library(readxl)
data_cattle_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_cattle_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_cattle_neg_2100.pdf", width = 4, height = 4)
treemap(data_cattle_neg,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Cattle", #指定面积大小的列
        vColor="Negatively-Impacted Cattle", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_cattle_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_cattle_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_cattle_neg_2100_Hunger.pdf", width = 4, height = 4)
treemap(data_cattle_neg,
        index=c("Hunger","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Cattle", #指定面积大小的列
        vColor="Negatively-Impacted Cattle", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()


data_cattle_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_cattle_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_cattle_neg_2100_FSI.pdf", width = 4, height = 4)
treemap(data_cattle_neg,
        index=c("FSI","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Cattle", #指定面积大小的列
        vColor="Negatively-Impacted Cattle", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()


data_cattle_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_cattle_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_cattle_neg_2100_Gender.pdf", width = 4, height = 4)
treemap(data_cattle_neg,
        index=c("Gender","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Cattle", #指定面积大小的列
        vColor="Negatively-Impacted Cattle", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

library(treemap)
library(readxl)
data_sheep_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_sheep_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_sheep_neg_2100.pdf", width = 4, height = 4)
treemap(data_sheep_neg,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Sheep", #指定面积大小的列
        vColor="Negatively-Impacted Sheep", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_sheep_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_sheep_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_sheep_neg_2100_Hunger.pdf", width = 4, height = 4)
treemap(data_sheep_neg,
        index=c("Hunger","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Sheep", #指定面积大小的列
        vColor="Negatively-Impacted Sheep", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()


data_sheep_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_sheep_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_sheep_neg_2100_FSI.pdf", width = 4, height = 4)
treemap(data_sheep_neg,
        index=c("FSI","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Sheep", #指定面积大小的列
        vColor="Negatively-Impacted Sheep", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()


data_sheep_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_sheep_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_sheep_neg_2100_Gender.pdf", width = 4, height = 4)
treemap(data_sheep_neg,
        index=c("Gender","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Sheep", #指定面积大小的列
        vColor="Negatively-Impacted Sheep", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

library(treemap)
library(readxl)
data_goats_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_goats_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_goats_neg_2100.pdf", width = 4, height = 4)
treemap(data_goats_neg,
        index=c("Development","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Goats", #指定面积大小的列
        vColor="Negatively-Impacted Goats", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()

data_goats_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_goats_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_goats_neg_2100_Hunger.pdf", width = 4, height = 4)
treemap(data_goats_neg,
        index=c("Hunger","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Goats", #指定面积大小的列
        vColor="Negatively-Impacted Goats", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()


data_goats_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_goats_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_goats_neg_2100_FSI.pdf", width = 4, height = 4)
treemap(data_goats_neg,
        index=c("FSI","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Goats", #指定面积大小的列
        vColor="Negatively-Impacted Goats", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()


data_goats_neg <- read_excel("/Users/lichaohui/Desktop/calculation/grazingniche/code/r_2100/treedata_goats_neg.xlsx")
pdf("/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treeplot_goats_neg_2100_Gender.pdf", width = 4, height = 4)
treemap(data_goats_neg,
        index=c("Gender","Country"), #指定多个分组的列，先按continent分组，再按country分组
        vSize="Negatively-Impacted Goats", #指定面积大小的列
        vColor="Negatively-Impacted Goats", #指定颜色深浅的列
        type="dens",
        title.legend = NA,  # Remove legend title
        position.legend = "none")  # Turn off legend)
dev.off()