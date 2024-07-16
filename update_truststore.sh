#!/bin/bash

GENIE_GENOMENEXUS_SITE="genie.genomenexus.org"
GENOMENEXUS_SERVERNAME="www.genomenexus.org"

truststore=$1
truststore_password=$2

if [ -z $truststore ] ; then
	echo "Missing truststore, exiting..."
	exit 1
fi

if [ -z $truststore_password ] ; then
	echo "Missing truststore password, exiting..."
	exit 1
fi

which openssl
if [ $? -ne 0 ] ; then
	echo "openssl needs to be installed. Exiting..."
	exit 1
fi

which keytool
if [ $? -ne 0 ] ; then
	echo "keytool needs to be installed. Exiting..."
	exit 1
fi

true | openssl s_client -connect $GENIE_GENOMENEXUS_SITE:443 -servername $GENOMENEXUS_SERVERNAME 2>/dev/null | openssl x509 > tmp.cert

keytool -delete -noprompt -alias $GENIE_GENOMENEXUS_SITE -keystore $truststore -storepass $truststore_password
keytool -import -noprompt -alias $GENIE_GENOMENEXUS_SITE -file tmp.cert -keystore $truststore -storepass $truststore_password
rm tmp.cert

echo "$truststore updated with newest cert from $GENOMENEXUS_SITE"
