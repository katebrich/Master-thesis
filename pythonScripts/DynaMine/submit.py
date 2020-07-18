#!/usr/bin/env python
# encoding: utf-8
"""
submit.py

Created by Elisa Cilia on 2014-08-20.
Copyright (c) 2014 Elisa Cilia. All rights reserved.
"""

import json
import urllib
#import urllib2
import time
import sys
import os

from . config import *


class DynaMineJSONInterface:

    def __init__(self, json_api_key):
        self._running = False
        self._json_api_key = json_api_key

    def submit_sequences(self, proteins, predictions_only = True):
        job = {'protocol': '1.0',
               'json_api_key': self._json_api_key,
               'sequences': proteins,
               'predictions_only': predictions_only,
        }
        return self._submit_job(job)

    def submit_uniprot_ids(self, uniprot_ids, predictions_only = True):
        job = {'protocol': '1.0',
               'json_api_key': self._json_api_key,
               'uniprot_ids': uniprot_ids,
               'predictions_only': predictions_only,
        }
        return self._submit_job(job)

    def _submit_job(self, job):
        # response gives us a job_id
        response = self._dynamine_request(job)
        if response['status'] == 'error':
            sys.stderr.write(response['message'] + '\n')
            sys.exit(1)
        job_id = response['job_id']
        self._print_progress(response['status'], False)
        # this makes the call blocking until the results are ready
        while response['status'] != 'completed':
            time.sleep(1.5)
            response = self._poll_results(job_id)
            if response['status'] == 'error':
                sys.stderr.write(response['message'] + '\n')
                sys.exit(1)
            self._print_progress(response['status'], True)
        return response['results']

    def _dynamine_request(self, request):
        # this is the post request with the json encoding of the job
        data = urllib.parse.urlencode({'batch': json.dumps(request)}).encode("utf-8")
        url = os.path.join(configuration_file['server_name'],'batch_request')
        response = ''
        try:
            #req = urllib2.Request(url, data)
            #http_response = urllib2.urlopen(req)
            http_response = urllib.request.urlopen(url, data)
            response = json.loads(http_response.read().decode('utf-8'))
        except Exception as e:
            response = {'status': 'error', 'message': 'Unable to communicate \
with the server. Please verify that %s is up and running, otherwise concact \
the server administrator (ERROR: %s).'
% (configuration_file['server_name'], str(e))}
        return response

    def _poll_results(self, job_id):
        request = {'protocol': '1.0',
               'json_api_key': self._json_api_key,
               'job_id': job_id
        }
        return self._dynamine_request(request)

    def _print_progress(self, status, flag):
        if status == 'queued':
            if flag:
                sys.stdout.write('Waiting to be processed.\n')
            else:
                sys.stdout.write('The request has been submitted.\n')
            self._running = False
        elif status == 'running':
            if not self._running:
                self._running = True
                sys.stdout.write('The request is being processed...')
            else:
                sys.stdout.write('.')
        elif status == 'completed':
            if not self._running:
                sys.stdout.write('The request is being processed...')
            sys.stdout.write('done.\n')
        elif status == 'error':
            sys.stderr.write('\n\nError while processing the request:\n\n')
        sys.stdout.flush()


if __name__ == '__main__':
    json_api_key = 'd03e71ca849d3bc5277fd68e7a0f45124df5be84bf348db65218982f'
    dynamine = DynaMineJSONInterface(json_api_key)

    proteins = {'sp|P04637|P53_HUMAN':
'''MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD''',
                'sp|Q92793|CBP_HUMAN':
'''MAENLLDGPPNPKRAKLSSPGFSANDSTDFGSLFDLENDLPDELIPNGGELGLLNSGNLV
PDAASKHKQLSELLRGGSGSSINPGIGNVSASSPVQQGLGGQAQGQPNSANMASLSAMGK
SPLSQGDSSAPSLPKQAASTSGPTPAASQALNPQAQKQVGLATSSPATSQTGPGICMNAN
FNQTHPGLLNSNSGHSLINQASQGQAQVMNGSLGAAGRGRGAGMPYPTPAMQGASSSVLA
ETLTQVSPQMTGHAGLNTAQAGGMAKMGITGNTSPFGQPFSQAGGQPMGATGVNPQLASK
QSMVNSLPTFPTDIKNTSVTNVPNMSQMQTSVGIVPTQAIATGPTADPEKRKLIQQQLVL
LLHAHKCQRREQANGEVRACSLPHCRTMKNVLNHMTHCQAGKACQVAHCASSRQIISHWK
NCTRHDCPVCLPLKNASDKRNQQTILGSPASGIQNTIGSVGTGQQNATSLSNPNPIDPSS
MQRAYAALGLPYMNQPQTQLQPQVPGQQPAQPQTHQQMRTLNPLGNNPMNIPAGGITTDQ
QPPNLISESALPTSLGATNPLMNDGSNSGNIGTLSTIPTAAPPSSTGVRKGWHEHVTQDL
RSHLVHKLVQAIFPTPDPAALKDRRMENLVAYAKKVEGDMYESANSRDEYYHLLAEKIYK
IQKELEEKRRSRLHKQGILGNQPALPAPGAQPPVIPQAQPVRPPNGPLSLPVNRMQVSQG
MNSFNPMSLGNVQLPQAPMGPRAASPMNHSVQMNSMGSVPGMAISPSRMPQPPNMMGAHT
NNMMAQAPAQSQFLPQNQFPSSSGAMSVGMGQPPAQTGVSQGQVPGAALPNPLNMLGPQA
SQLPCPPVTQSPLHPTPPPASTAAGMPSLQHTTPPGMTPPQPAAPTQPSTPVSSSGQTPT
PTPGSVPSATQTQSTPTVQAAAQAQVTPQPQTPVQPPSVATPQSSQQQPTPVHAQPPGTP
LSQAAASIDNRVPTPSSVASAETNSQQPGPDVPVLEMKTETQAEDTEPDPGESKGEPRSE
MMEEDLQGASQVKEETDIAEQKSEPMEVDEKKPEVKVEVKEEEESSSNGTASQSTSPSQP
RKKIFKPEELRQALMPTLEALYRQDPESLPFRQPVDPQLLGIPDYFDIVKNPMDLSTIKR
KLDTGQYQEPWQYVDDVWLMFNNAWLYNRKTSRVYKFCSKLAEVFEQEIDPVMQSLGYCC
GRKYEFSPQTLCCYGKQLCTIPRDAAYYSYQNRYHFCEKCFTEIQGENVTLGDDPSQPQT
TISKDQFEKKKNDTLDPEPFVDCKECGRKMHQICVLHYDIIWPSGFVCDNCLKKTGRPRK
ENKFSAKRLQTTRLGNHLEDRVNKFLRRQNHPEAGEVFVRVVASSDKTVEVKPGMKSRFV
DSGEMSESFPYRTKALFAFEEIDGVDVCFFGMHVQEYGSDCPPPNTRRVYISYLDSIHFF
RPRCLRTAVYHEILIGYLEYVKKLGYVTGHIWACPPSEGDDYIFHCHPPDQKIPKPKRLQ
EWYKKMLDKAFAERIIHDYKDIFKQATEDRLTSAKELPYFEGDFWPNVLEESIKELEQEE
EERKKEESTAASETTEGSQGDSKNAKKKNNKKTNKNKSSISRANKKKPSMPNVSNDLSQK
LYATMEKHKEVFFVIHLHAGPVINTLPPIVDPDPLLSCDLMDGRDAFLTLARDKHWEFSS
LRRSKWSTLCMLVELHTQGQDRFVYTCNECKHHVETRWHCTVCEDYDLCINCYNTKSHAH
KMVKWGLGLDDEGSSQGEPQSKSPQESRRLSIQRCIQSLVHACQCRNANCSLPSCQKMKR
VVQHTKGCKRKTNGGCPVCKQLIALCCYHAKHCQENKCPVPFCLNIKHKLRQQQIQHRLQ
QAQLMRRRMATMNTRNVPQQSLPSPTSAPPGTPTQQPSTPQTPQPPAQPQPSPVSMSPAG
FPSVARTQPPTTVSTGKPTSQVPAPPPPAQPPPAAVEAARQIEREAQQQQHLYRVNINNS
MPPGRTGMGTPGSQMAPVSLNVPRPNQVSGPVMPSMPPGQWQQAPLPQQQPMPGLPRPVI
SMQAQAAVAGPRMPSVQPPRSISPSALQDLLRTLKSPSSPQQQQQVLNILKSNPQLMAAF
IKQRTAKYVANQPGMQPQPGLQSQPGMQPQPGMHQQPSLQNLNAMQAGVPRPGVPPQQQA
MGGLNPQGQALNIMNPGHNPNMASMNPQYREMLRRQLLQQQQQQQQQQQQQQQQQQGSAG
MAGGMAGHGQFQQPQGPGGYPPAMQQQQRMQQHLPLQGSSMGQMAAQMGQLGQMGQPGLG
ADSTPNIQQALQQRILQQQQMKQQIGSPGQPNPMSPQQHMLSGQPQASHLPGQQIATSLS
NQVRSPAPVQSPRPQSQPPHSSPSPRIQPQPSPHHVSPQTGSPHPGLAVTMASSIDQGHL
GNPEQSAMLPQLNTPSRSALSSELSLVGDTTGDTLEKFVEGL'''
    }
    results = dynamine.submit_sequences(proteins, predictions_only = False)
    print(results)

    uniprot_ids = ['P04637', 'P04638']
    results = dynamine.submit_uniprot_ids(uniprot_ids, predictions_only = True)
    print(results)


