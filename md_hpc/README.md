
# About

Esse `Dockerfile` consegue executar as etapas de MD e análise utilizando `autoGromacs` em ambiente HPC (NVidia).
A construção da imagem utiliza os campos de força de GROMACS contidos na pasta `FF`.

O `Makefile` contém operações comuns da manutenção da imagem `docker`.

Os comandos exatos a serem executados dentro da imagem estão disponíveis no `autogromacs --help` ou no repositório daquele software.
