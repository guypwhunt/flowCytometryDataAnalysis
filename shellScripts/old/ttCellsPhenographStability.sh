#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH -p cpu
#SBATCH --mem=1000G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

source /etc/profile.d/modules.sh
module load bzip2/1.0.8-gcc-9.4.0 zlib/1.2.11-gcc-9.4.0 boost/1.79.0-gcc-9.4.0-pic-no-python icu4c/67.1-gcc-9.4.0 linux-pam/1.5.2-gcc-9.4.0 openjdk/11.0.12_7-gcc-9.4.0 font-util/1.3.2-gcc-9.4.0 libiconv/1.16-gcc-9.4.0 xz/5.2.5-gcc-9.4.0 libxml2/2.9.12-gcc-9.4.0 util-linux-uuid/2.36.2-gcc-9.4.0 fontconfig/2.13.94-gcc-9.4.0 ncurses/6.2-gcc-9.4.0 tar/1.34-gcc-9.4.0 gettext/0.21-gcc-9.4.0 libffi/3.3-gcc-9.4.0 pcre/8.44-gcc-9.4.0 berkeley-db/18.1.40-gcc-9.4.0 readline/8.1-gcc-9.4.0 gdbm/1.19-gcc-9.4.0 libmd/1.0.3-gcc-9.4.0 libbsd/0.11.3-gcc-9.4.0 expat/2.4.1-gcc-9.4.0 openssl/1.1.1l-gcc-9.4.0 sqlite/3.36.0-gcc-9.4.0 python/3.8.12-gcc-9.4.0 glib/2.70.0-gcc-9.4.0-python-3.8.12 inputproto/2.3.2-gcc-9.4.0 kbproto/1.0.7-gcc-9.4.0 libpthread-stubs/0.4-gcc-9.4.0 libxau/1.0.8-gcc-9.4.0 libxdmcp/1.1.2-gcc-9.4.0 xcb-proto/1.14.1-gcc-9.4.0 xextproto/7.3.0-gcc-9.4.0 xproto/7.0.31-gcc-9.4.0 xtrans/1.3.5-gcc-9.4.0 libx11/1.7.0-gcc-9.4.0 libxcb/1.14-gcc-9.4.0 libxext/1.3.3-gcc-9.4.0 renderproto/0.11.1-gcc-9.4.0 libxrender/0.9.10-gcc-9.4.0 cairo/1.16.0-gcc-9.4.0-python-3.8.12 freetype/2.11.0-gcc-9.4.0 libice/1.0.9-gcc-9.4.0 libsm/1.2.3-gcc-9.4.0 dbus/1.12.8-gcc-9.4.0-python-3.8.12 libxfixes/5.0.2-gcc-9.4.0 libxtst/1.2.2-gcc-9.4.0 recordproto/1.14.2-gcc-9.4.0 at-spi2-core/2.40.1-gcc-9.4.0-python-3.8.12 at-spi2-atk/2.38.0-gcc-9.4.0-python-3.8.12 atk/2.36.0-gcc-9.4.0-python-3.8.12 fixesproto/5.0-gcc-9.4.0 fribidi/1.0.5-gcc-9.4.0 gdk-pixbuf/2.42.6-gcc-9.4.0-python-3.8.12 gobject-introspection/1.56.1-gcc-9.4.0-python-3.8.12 glproto/1.4.17-gcc-9.4.0 mesa18/18.3.6-gcc-9.4.0 libepoxy/1.4.3-gcc-9.4.0 libcroco/0.6.13-gcc-9.4.0-python-3.8.12 graphite2/1.3.13-gcc-9.4.0 libxft/2.3.2-gcc-9.4.0 librsvg/2.51.0-gcc-9.4.0-python-3.8.12 libxi/1.7.6-gcc-9.4.0 util-macros/1.19.3-gcc-9.4.0 xkbdata/1.0.1-gcc-9.4.0 libxkbcommon/0.8.2-gcc-9.4.0 pango/1.42.0-gcc-9.4.0-python-3.8.12 shared-mime-info/1.9-gcc-9.4.0-python-3.8.12 libxrandr/1.5.0-gcc-9.4.0 randrproto/1.5.0-gcc-9.4.0 xrandr/1.5.0-gcc-9.4.0 gtkplus/3.24.29-gcc-9.4.0-python-3.8.12 krb5/1.19.3-gcc-9.4.0 lcms/2.9-gcc-9.4.0 libjpeg-turbo/2.1.0-gcc-9.4.0 libtiff/4.3.0-gcc-9.4.0 ghostscript/9.54.0-gcc-9.4.0-python-3.8.12 gmp/6.2.1-gcc-9.4.0 harfbuzz/2.6.8-gcc-9.4.0-python-3.8.12 libgd/2.2.4-gcc-9.4.0 libpaper/1.1.28-gcc-9.4.0 libpng/1.6.37-gcc-9.4.0 libxmu/1.1.2-gcc-9.4.0 libxpm/3.5.12-gcc-9.4.0 libxaw/1.0.13-gcc-9.4.0 libxt/1.1.5-gcc-9.4.0 mpfr/4.1.0-gcc-9.4.0 perl/5.34.0-gcc-9.4.0 pixman/0.40.0-gcc-9.4.0 poppler-data/0.4.9-gcc-9.4.0 poppler/0.79.0-gcc-9.4.0 teckit/2.5.9-gcc-9.4.0 zziplib/0.13.72-gcc-9.4.0-python-3.8.12 texlive/20210325-gcc-9.4.0-python-3.8.12 pandoc/2.14.0.3-gcc-9.4.0-python-3.8.12 curl/7.79.0-gcc-9.4.0 pcre2/10.36-gcc-9.4.0 scrnsaverproto/1.2.2-gcc-9.4.0 libxscrnsaver/1.2.2-gcc-9.4.0 tcl/8.6.11-gcc-9.4.0 tk/8.6.11-gcc-9.4.0 r/4.1.1-gcc-9.4.0-withx-rmath-standalone-python-3.8.12 postgresql/14.0-gcc-9.4.0 soci/4.0.2-gcc-9.4.0 yaml-cpp/0.6.3-gcc-9.4.0 yarn/1.22.4-gcc-9.4.0 rstudio_kcl/v2022.07.1_554-gcc-9.4.0-r4.1.1-python-3.8.12 

iteration=${1}
addition=1
limit=8

if [ "$iteration" -lt "$limit" ]; then
  echo "Iteration is less than limit"
  echo $iteration
  Rscript R/tCells/21_phenograph_bootstrapping.R
  sbatch shellScripts/ttCellsPhenographStability.sh $((iteration + addition));
fi