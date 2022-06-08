

CNVscore: A framework for the prioritization of CNVs with uncertainty estimates in rare disease patients

-----

Welcome to the official Github repository of **CNVscore**.


The CNVscore manuscript is available on:

TBA


# Table of contents

  - [Overview](#Overview)
  - [API query example](#api-query-example)
  - [Response example](#response-example)
  - [Availability](#Availability)
  - [Docker installation](#docker-installation)
  - [Authors and contact](#authors-and-contact)
  - [License](#License)
  - [Disclaimer](#Disclaimer)
  - [Reference](#Reference)
  - [News](#News)

## Overview

[insert]

<p align="center">

<img  src="https://github.com/RausellLab/CNVscore/blob/main/doc/CNVscore_overview.svg">

</p>

## API query example

``` bash
TBA
```

## Response example

``` bash
[
  {
    "chrom": "17",
    "start": "61,158,866",
    "end": "61,353,248",
    "variant_class": "deletion",
    "cnvscore": 0.9132,
    "uncertainty_level": 1,
    "rules": "min_expression<=39 & max_cadd>34.5 (risk:0.87 - support:3910);min_expression<=33.5 & max_cadd>37.5 (risk:0.88 - support:3445);loeuf>0.8085 & cpg_density>1.5 (risk:0.81 - support:789);enhancer<=0.000282170504594472 & max_cadd>40.5 (risk:0.87 - support:3544);pli>82.5 & max_cadd>25.7 (risk:0.91 - support:1385)"
  }
]
```

## Availability

CNVscore can be deployed as a private API service through a Docker image without external dependencies. Instructions to locally deploy the API service are provided in the next section.

In addition, CNVscore can be queried and interrogated at <http://cnvxplorer.com>. 



## Docker installation

``` bash

# Note: the first session after the deployment is slower since the application loads all the data required

git clone https://github.com/RausellLab/CNVscore.git

cd CNVscore 

docker build -t cnvscore . -f api/Dockerfile . # The tag "cnvscore" is optional

docker run -d -p 3838:3838 cnvscore # -p (specify port) -d (detached mode)

# The port 3838 is optional. Please make sure you set a port not blocked by firewalls.
# If you change the port number (3838) by any other, make sure to set it in the Dockerfile (EXPOSE instruction)
```

## Authors and contact

CNVscore has been developed by Francisco Requena and Antonio Rausell,
at the [Clinical Bioinformatics
Laboratory](https://www.institutimagine.org/en/antonio-rausell-161) of
the [Imagine Institute](https://www.institutimagine.org/en/) in Paris,
France.

Please address comments and questions about CNVscore to: \*
**Francisco Requena** -
[francisco.requena@institutimagine.org](francisco.requena@institutimagine.org)
\* **Antonio Rausell** -
[antonio.rausell@institutimagine.org](antonio.rausell@institutimagine.org)

## License

This project is licensed under the GNU General Public License 3 - see
the [LICENSE](LICENSE) file for details

See the License for the specific language governing permissions and
limitations under the License.

Copyright 2021 Clinical BioInformatics Laboratory - Institut Imagine

## Disclaimer

CNVscore or any document available from this server are distributed on
an “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
express, implied, or statutory, including, but not limited to, any
implied warranties of merchantability, fitness for a particular purpose
and freedom from infringement, or that CNVscore or any documents
available from this server will be error-free.

In no event will the Imagine Institute, the Clinical Bioinformatics lab,
or any of its members be liable for any damages, including but not
limited to direct, indirect, special, or consequential damages, arising
out of, resulting from, or in any way connected with the use of
CNVscore or documents available from it.

## Reference

TBA

## News

You may follow us in Twitter for regular news and updates:
<https://twitter.com/AntonioRausell>
