{{/* Generate a template that can be used for both assigned and unassigned xperiments */}}
# go templating language
{{- define "pipeline.pod-template" -}}
    metadata:
      name: "{{ .Release.Name }}-server"
      labels:
        type: 'pipeline'
        sandboxId: "{{ .Values.sandboxId }}"
    spec:
      restartPolicy: Always
      serviceAccountName: 'deployment-runner'
      containers:
      - name: "{{ .Release.Name }}"
        image: "{{ .Values.pipelineRunner.image }}"
        env:
          - name: CLUSTER_ENV
            value: "{{ .Values.clusterEnv }}"
          - name: SANDBOX_ID
            value: "{{ .Values.sandboxId }}"
          - name: AWS_ACCOUNT_ID
            value: "{{ .Values.myAccount.accountId }}"
        volumeMounts:
        - name: podinfo
          mountPath: /etc/podinfo
        resources:
          requests:
            memory: "{{ .Values.memoryRequest }}"
      - name: datadog-agent
        image: datadog/agent
        env:
        - name: DD_API_KEY
          value: "{{ .Values.myAccount.datadogApiKey }}"
        - name: DD_SITE
          value: "datadoghq.eu"
        - name: DD_EKS_FARGATE
          value: "true"
        - name: DD_CLUSTER_NAME
          value: "biomage-{{ .Values.clusterEnv }}"
        # Disable log collection by DD agent
        # because we use FluentBit for Fargate log collection
        - name: DD_LOGS_ENABLED
          value: "false"
        - name: DD_CONTAINER_EXCLUDE
          value: "name:.*"
        - name: DD_KUBERNETES_POD_LABELS_AS_TAGS
          value: '{"*": "%%label%%"}'
        - name: DD_KUBERNETES_KUBELET_NODENAME
          valueFrom:
            fieldRef:
              apiVersion: v1
              fieldPath: spec.nodeName
      volumes:
      - name: podinfo
        downwardAPI:
          items:
            - path: "labels"
              fieldRef:
                fieldPath: metadata.labels
{{- end -}}
