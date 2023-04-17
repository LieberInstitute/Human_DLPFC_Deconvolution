library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("readxl")
library("ggrepel")
library("jaffelab")
library("plotly")

## prep dirs ##
plot_dir <- here("plots", "02_quality_control", "02_bulk_qc_plotly")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load New Data ####
# load(here("processed-data","01_SPEAQeasy","round1_2022-07-06","rse_gene.Rdata"), verbose = TRUE)
load(here("processed-data", "01_SPEAQeasy", "round2_v40_2022-07-06", "rse", "rse_gene.Rdata"), verbose = TRUE)
pd_new <- as.data.frame(colData(rse_gene))

## Subset key QC metrics
qc_variables <- c("numReads", "numMapped", "numUnmapped", "overallMapRate", "concordMapRate", "totalMapped", "mitoMapped", "mitoRate", "rRNA_rate", "totalAssignedGene")

pd_new_qc_long <- pd_new |>
    select(SAMPLE_ID, seq_set, library_type, library_prep, all_of(qc_variables)) |>
    pivot_longer(!c(SAMPLE_ID, seq_set, library_type, library_prep), names_to = "qc_var")

head(pd_new_qc_long)

#### plotly boxplots ####
pd_new_qc_key <- pd_new_qc_long |>
    highlight_key(~SAMPLE_ID)

qc_boxplots <- pd_new_qc_key |>
    ggplot(aes(x = seq_set, y = value)) +
    geom_boxplot(aes(fill = `library_type`), outlier.shape = NA) +
    geom_jitter() +
    facet_wrap(~qc_var, scales = "free_y", ncol = 5) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17)
    ) +
    labs(title = "Human DLPFC Deconvolution", subtitle = "Bulk Round 1 + 2")

ggsave(qc_boxplots, filename = here(plot_dir, "qc_boxplots.png"), width = 20, height = 12)

qc_boxplots_plotly <- ggplotly(qc_boxplots, tooltip = c("colour", "text"))

htmlwidgets::saveWidget(
    highlight(qc_boxplots_plotly,
        on = "plotly_click",
        off = "plotly_doubleclick",
        selectize = TRUE,
        dynamic = TRUE,
        persistent = FALSE,
    ),
    selfcontained = FALSE,
    file = here(plot_dir, "qc_boxplots.html")
)


#### QC cutoffs  ####
qc_cutoff_tb <- tibble(
    metric = c("overallMapRate", "totalAssignedGene", "numReads", "rRNA_rate"),
    cutoff = c(0.5, 0.3, 10^7.25, 1e-3)
)

## test static plot
pd_new_qc_long |>
    group_by(library_type) |>
    group_map(~ {
        qc_boxplots_lt <- ggplot(data = .x, aes(x = seq_set, y = value)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter() +
            # geom_hline(data = qc_cutoff_tab, y_intercept = cutoff)+
            facet_wrap(~qc_var, scales = "free_y", ncol = 5) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                text = element_text(size = 17)
            ) +
            labs(title = paste("DLPFC Bulk: ", .y))

        ggsave(qc_boxplots_lt, filename = here(plot_dir, paste0("qc_boxplots_", .y, ".png")), width = 20, height = 12)
    })

## Plotly by libary_type ###
pd_new_qc_long |>
    group_by(library_type) |>
    group_map(~ {
        qc_key_lt <- .x |>
            highlight_key(~SAMPLE_ID)

        qc_boxplots_lt <- ggplot(data = qc_key_lt, aes(x = seq_set, y = value)) +
            geom_violin() +
            geom_jitter() +
            facet_wrap(~qc_var, scales = "free_y", ncol = 5) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                text = element_text(size = 17)
            ) +
            labs(title = paste("DLPFC Bulk: ", .y))

        # ggsave(qc_boxplots_lt, filename = here(plot_dir, paste0("qc_boxplots_",.y,".png")), width = 20, height = 12)

        qc_boxplots_lt_plotly <- ggplotly(qc_boxplots_lt, tooltip = c("colour", "text"))

        htmlwidgets::saveWidget(
            highlight(qc_boxplots_lt_plotly,
                on = "plotly_click",
                off = "plotly_doubleclick",
                selectize = TRUE,
                dynamic = TRUE,
                persistent = FALSE,
            ),
            selfcontained = FALSE,
            file = here(plot_dir, paste0("qc_boxplots_", .y, ".html"))
        )
    })



pd_new |> filter(totalMapped < 5e7)





# sgejobs::job_single('01_prelim_bulk_qc_check', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 01_prelim_bulk_qc_check.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
