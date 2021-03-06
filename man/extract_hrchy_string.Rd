% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_hrchy_string.R
\name{extract_hrchy_string}
\alias{extract_hrchy_string}
\title{Extract hierarchy from a string (based on prior knowledge)}
\usage{
extract_hrchy_string(input_string, special_character = "@", plot = T)
}
\arguments{
\item{input_string}{An input string with nested labeled parentheses that characterizes hierarchy of annotated cell types. Avoid using comma (','), parentheses ('()') and special character (e.g.'@') in celltype names.}

\item{special_character}{A character denotes non-leaf node (including root node and internal nodes).}

\item{plot}{Plot tree structure or not.}
}
\value{
A list of hierarchy information:
        \itemize{
        \item edges: a data frame that documents that parent-children relationship.
        \item layout: a data frame that contains id, xy coordiates, label and a leaf indicator for each node.
        \item immediate_children: a data frame that captures immediate children for each node.
        \item leaves_info: a data frame that summarizes leaf children for each node.
        }
}
\description{
Parse an input string by extracting parent-children relationship.
}
\examples{
input_string = '@T(@CD4+(CD4_c1,CD4_c2),@CD8+(CD8_c1,CD8_c2))'
extract_hrchy_string(input_string, special_character = '@')
}
\author{
Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
}
