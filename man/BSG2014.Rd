\name{BSG2014}
\alias{BSG2014}
\docType{data}
\title{
Business simulation game data
}
\description{
The data were collected from 354 senior business administration students during
a business simulation game at a Western European University.

The game was played for a total of 12 rounds (i.e., two separate games of 6
rounds) as part of the capstone strategy class.  Students were randomly
assigned to teams of four, and surveyed in three waves: prior to the first
game, in between the two games, and after the second game (with different
variables being surveyed in the different waves).

The 354 students formed 92 teams, and the responses of individual students were
aggregated to the team level.  Leaving out teams with less than 50 percent
response rate yields \eqn{n = 89} teams.  Only a small subset of the collected
variables are included here.
}
\usage{data("BSG2014")}
\format{
  A data frame with 89 observations on the following 13 variables.
  \describe{

    \item{\code{ProcessConflict}}{Process conflict was operationalized with the
    three-item scale of Jehn (1995) and responses were aggregated.  Process
    conflict was measured in the second survey (between the two games).}

    \item{\code{SharedExperience}}{Teams were randomly formed, no prior shared
    group experience is expected and shared group experience and training is
    developed during the first game for the second game.  The team's score on
    the first game, which is determined by objective performance measures, is
    used as a proxy for the level of shared group experience and training.}

    \item{\code{TaskConflict}}{Task conflict was operationalized with the
    intra-group conflict scale of Jehn (1995).  Five items on the presence of
    conflict were rated on a 5-point Likert scale (1 = none, 5 = a lot) and
    aggregated.  Task conflict was measured in the second survey (between the
    two games).}

    \item{\code{TeamCommitment}}{Team commitment was measured by four items
    based on Mowday, Steers & Porter (1979) and responses were aggregated.
    Team commitment was measured in the third survey (after the second game).}

    \item{\code{TeamPerformance}}{Team performance in the second game was
    measured subjectively by the team members’ perceptions of the team’s
    functioning, using aggregated responses to Hackman’s (1986) Likert scale
    items.  Subjective team performance was measured in the third survey (after
    the second game).}

    \item{\code{TMS}}{Transactive memory systems (TMS) are defined as shared
    systems that people in relationships develop for encoding, storing, and
    retrieving information about different substantive domains.  TMS was
    operationalized with Lewis’ (2003) 15-item scale that measures the three
    sub-dimensions of TMS (credibility, specialization and coordination).  Team
    members responded on a 5-point scale (1 = strongly disagree, 5 = strongly
    agree).  Following Lewis (2003), the three sub dimensions were aggregated
    to form the TMS construct.  TMS was measured in the second survey (between
    the two games).}

    \item{\code{ValueDiversity}}{Value diversity was operationalized with the
    short version of Schwartz’s Value Survey (SVS) to measure team members’
    individual values (Lindeman & Verkasalo, 2005).  The responses were
    aggregated with the average of the coefficient of variations of each value
    dimension among team members.  Value diversity was measured in the first
    survey (before the first game).}

    \item{\code{ProceduralJustice}}{}

    \item{\code{InteractionalJustice}}{}

    \item{\code{SharedLeadership}}{}

    \item{\code{AgeDiversity}}{}

    \item{\code{GenderDiversity}}{}

    \item{\code{TeamScore}}{The team's score on the second game, which is
    determined by objective performance measures.}

  }
}
\references{
Hackman, J.R. (1986) The Psychology of Self-Management in Organizations.  In
Pallack, M.S and Perloff, R.O. (Eds.), \emph{Psychology and Work: Productivity,
Change, and Employment}, 89--136.  Washington, DC: American Psychological
Association.

Jehn, K.A. (1995) A Multi-Method Examination of the Benefits and Detriments of
Intra-Group Conflict.  \emph{Administrative Science Quarterly}, \bold{40}(2),
256--285.

Lewis, K. (2003) Measuring Transactive Memory Systems in the Field: Scale
Development and Validation.  \emph{Journal of Applied Psychology},
\bold{88}(4), 587--604.

Lindeman, M. and Verkasalo, M. (2005) Measuring Values With the Short
Schwartz's Value Survey.  \emph{Journal of Personality Assessment},
\bold{85}(2), 170--178.

Mowday, R.T., Steers, R.M. and Porter, L.W. (1979) The Measurement of
Organizational Commitment.  \emph{Journal of Vocational Behavior},
\bold{14}(2), 224--47.
}
\examples{
data("BSG2014")
summary(BSG2014)

# scatterplot matrix for the variables included in the
# illustrative mediation analysis
x <- "ValueDiversity"
y <- "TeamCommitment"
m <- "TaskConflict"
plot(BSG2014[, c(x, y, m)], pch = 21, bg = "black")
}
\keyword{datasets}
