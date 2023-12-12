make_gene_panel <- function(chr, range) {
  require(gggenes)

  regions <- GRanges(seqnames = chr, ranges = IRanges(start = range[1], width = diff(range)))
  ann <- import('sequencing/GCF_000633615.1_Guppy_female_1.0_MT_genomic.gff.gz', which = regions) %>%
    as.data.frame()

  genes <- ann %>% filter(type == 'gene')
  exons <- ann %>% filter(type == 'exon') %>%
    dplyr::select(substart = start, subend = end, gene) %>%
    left_join(dplyr::select(genes, start, end, gene, strand), 'gene')

  genes$gene <- ifelse(genes$gene == 'LOC103476393', 'texim', genes$gene)
  exons$gene <- ifelse(exons$gene == 'LOC103476393', 'texim', exons$gene)

  if (nrow(genes) == 0) {
    p <- ggplot() + theme_void()
  } else {
    p <- ggplot() +
      geom_gene_arrow(
        aes(xmin = start, xmax = end, y = gene, forward = strand == '+'),
        genes, fill = 'grey90'
      ) +
      geom_subgene_arrow(
        aes(xmin = start, xmax = end, xsubmin = substart, xsubmax = subend, y = gene, forward = strand == '+'),
        exons, fill = 'grey70', color = NA
      ) +
      geom_text(
        aes(x = (start + end) / 2, y = gene, label = gene),
        genes, size = 3
      ) +
      scale_y_discrete(expand = expansion(mult = c(0.2, 0.5))) +
      scale_x_continuous(
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      ) +
      coord_cartesian(xlim = range) +
      theme_void()
  }
  return(p)
}
