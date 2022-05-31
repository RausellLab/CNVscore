

CNVscore: A framework for the prioritization of CNVs with uncertainty estimates in rare disease patients

-----

Welcome to the official Github repository of **CNVscore**.


The CNVscore manuscript is available on:

[insert]


# Table of contents

  - [Overview](#Overview)
  - [Availability](#Availability)
  - [Docker installation](#docker-installation)
  - [Authors and contact](#authors-and-contact)
  - [License](#License)
  - [Disclaimer](#Disclaimer)
  - [References](#References)
  - [News](#News)

## Overview

[insert]

<p align="center">

<img  src="https://github.com/RausellLab/CNVxplorer/blob/master/doc/Overview.svg">

</p>

## Availability

CNVscore can be deployed as a private API service through a Docker image without external dependencies. Instructions to locally deploy the API service are provided in the next section.

In addition, CNVscore can be queried and interrogated at <http://cnvxplorer.com>. 

## Example API query

### Individual CNV

``` bash
```

### Multiple CNVs

``` bash
```


## Docker installation

``` bash

# Note: the first session after the deployment is slower since the application loads all the data required

git clone https://github.com/RausellLab/CNVxplorer.git

mv CNVscore/Dockerfile .

docker build -t cnvscore . # The tag "cnvscore" is optional

docker run -d -p 3838:3838 cnvscore # -p (specify port) -d (detached mode)

# The port 3838 is optional. Please make sure you set a port not blocked by firewalls.
# If you change the port number (3838) by any other, make sure to set it in the Dockerfile 
# (EXPOSE instruction)
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

## References

[insert]

## News

ou may follow us in Twitter for regular news and updates:
<https://twitter.com/AntonioRausell>

# CNVscore
Machine-learning model for the prioritization of CNVs with uncertainty estimates in rare disease patients

