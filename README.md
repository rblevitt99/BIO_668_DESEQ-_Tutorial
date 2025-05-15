# BIO_668_DESEQ-_Tutorial
Here is a DESEQ2 tutorial I followed in my BIO 668 Course (bioinformatics)


```{r}
library(DESeq2) 
library(tidyverse)
library(ggplot2)
count=read.csv("rna_counts_data.csv")
count=subset(count,select= -c(Length))
View(count)

map=read.csv("rna_map_update copy.csv")
meta=c("Sample.Name","Condition")
metadata=map[meta]
View(metadata)
View(count)
count <- count %>% remove_rownames %>% column_to_rownames(var="Geneid")
View(count)
View(map)
```

```{r}
SampleName = c(colnames(count)) # obtaining sample names from the count matric
condition = c("Mutation","Mutation","Mutation","Mutation","Mutation","Mutation","Mutation","Mutation","WT","WT","Mutation","Mutation","Mutation","Mutation","Mutation","Mutation","Mutation","Mutation","Mutation","WT","WT","WT") # specifying the conditions
meta_data = data.frame(SampleName, condition) # making the metadata frame
meta_data
meta_data <- meta_data %>% remove_rownames %>%  column_to_rownames(var="SampleName") # making the sample name row ID
meta_data

```


```{r}
all(colnames(count) == rownames(meta_data)) # double checkong if the names in count match the names in the meta data




```



```{r}
dds = DESeqDataSetFromMatrix(countData = count,
                             colData = meta_data,
                             design = ~condition) # constructs the deseq data object. It will stoire stuff like read counts, metadata, and the condition parameters
dds
```
```{r}
# filter any counts less than 10
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]
dds
```

```{r}
# set the reference to be the control (WT) as R will automatically go in alphabetical order.
dds$condition <- relevel(dds$condition, ref ='WT')
```


```{r}
# this code will obtain the normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
```

```{r}
dds <- DESeq(dds) # running the differential expression
res <- results(dds) # setting res (results) to the results ouptput from the differential expression
res
summary(res)
sum(res$padj < 0.01, na.rm=TRUE) # finding the number of DE genes that have a p value < 0.01
```
```{r}
res05 <- results(dds, alpha=0.05) # same thing as above just with a p value < 0.05.
summary(res05)
sum(res$padj < 0.05, na.rm=TRUE)
```


```{r}
data <- data.frame(res)
head(data)
View(data)
```

```{r}
rld <- rlog(dds)
plotPCA(rld)
```


```{r}
data <- data %>%
  mutate(
  Expression = case_when(log2FoldChange >= log(1) & padj <= 0.05 ~ "Up-regulated",
  log2FoldChange <= -log(1) & padj <= 0.05 ~ "Down-regulated",
  TRUE ~ "Unchanged")
  )
head(data)
```


```{r}
top <- 10
# we are getting the top 10 up and down regulated genes by filtering the column Up-regulated and Down-regulated and sorting by the adjusted p-value. 
top_genes <- bind_rows(
  data %>%
  filter(Expression == 'Up-regulated') %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(top),
  data %>%
  filter(Expression == 'Down-regulated') %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(top)
  )
# create a datframe just holding the top 10 genes
Top_Hits = head(arrange(data,pvalue),10)
Top_Hits
View(Top_Hits)
```

```{r}
data$label = if_else(rownames(data) %in% rownames(Top_Hits), rownames(data), "")
# basic plot
p1 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
geom_point( size = 2/5) +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"P-Value")) +
xlim(-4.5, 4.5)
p1
```

```{r}
# basic plot with line + red for p < 0.05
p2 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
geom_point(aes(color = Expression), size = 2/5) +
#geom_hline(yintercept= -log10(0.05), linetype="dashed", linewidth = .4) +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"P-Value")) +
scale_color_manual(values = c("firebrick3", "black", "firebrick3")) +
xlim(-4.5, 4.5) +
theme(legend.position = "none")
p2
```


```{r}
library(ggrepel)
p3 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
geom_point(aes(color = Expression), size = 2/5) +
# geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"P-Value")) +
scale_color_manual(values = c("firebrick3", "black", "firebrick3")) +
xlim(-4.5, 4.5) +
theme(legend.position = "none") +
geom_text_repel(aes(label = label), size = 2.5)
p3
```


```{r}
# plot with up/down
p4 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
geom_point(aes(color = Expression), size = 2/5) +
#geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"P-Value")) +
scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
xlim(-4.5, 4.5)
p4
```

```{r}
p5 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
geom_point(aes(color = Expression), size = 2/5) +
# geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"P-Value")) +
scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
xlim(-4.5, 4.5) +
geom_text_repel(aes(label = label), size = 2.5)
p5
```
